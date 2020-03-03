###This is a file operation that converts summarised SNP data into a string
###of 1s and 0s that can be analysed by the HMM Viterbi algorithm.

###It takes in a file assigning a type to SNPs, discards irrelevant ones, and
###classifies relevant ones according to a chosen scheme.
###It writes out the SNPs as a single string and their positions as a key file.

###This is based on numpy arrays and should be way faster than fileoperation_classify_archaic_annotation_as_string

###GSJ 03/10/2017

import numpy as np
import copy
import gzip

def fileoperation_classify_archaic_annotation_as_string_array(infile, outfile_base, classification_scheme, exclusion_scheme = None, check_conflicts = True, freq_columns = []):
    """
    Takes in a file assigning a type to SNPs, discards irrelevant ones, and classifies relevant ones according to instructions in classification_scheme.
    classification_scheme is a list of conditions, with [[conditions_0], [conditions_1], ..., [conditions_N]] indicating how each SNP should be classified.
    conditions_0 is a nested list of alternative AND conditions:
    [[[column, value] AND [column, value]] OR [[column, value] AND [column, value]] etc.
    this should be slow but effective.
    The program will check_conflicts by default, but will be faster if this is disabled.
    There is also the option of suggesting an exclusion scheme, which is applied to the snp first and may speed things up by filtering out snps that are not required.
    V2 - numpy array based. Should be faster.
    V3 - specified columns can be format [column, min_value, max_value] to allow frequency filtering
    """
    open_f = gzip.open if infile[-3:] == '.gz' else open
    snp_pos = []
    snp_types = []
    conditions_list = []
    #The first step is to convert the classification system to a series of functions acting on numpy arrays.
    if len(classification_scheme) > 42:
        raise RuntimeError("Only implemented to work with types 0 to 42, with types > 10 being other characters.")
    snp_data = []
    headers = None
    with open_f(infile, 'rb') as f_in:
        for line in f_in:
            split_line = line[0:-1].split()
            if headers is None:
                headers = split_line
                pass
            else:
                snp_data.append(split_line)
    snp_data = np.array(snp_data)
    #Now, exclude some
    exclude = np.zeros(len(snp_data), dtype = bool)
    if exclusion_scheme is None:
        pass
    else:
        for exclusion_option in exclusion_scheme:
            exclusion_equals = []
            exclusion_freq = []
            for req in exclusion_option:
                if req[0] in freq_columns:
                    exclusion_freq.append(req)
                else:
                    exclusion_equals.append(req)
            equals_mask = np.product([snp_data[::,req[0]] == str(req[1]) for req in exclusion_equals], axis = 0, dtype = bool)
            if len(exclusion_freq) > 0:
                product_mask = np.product([np.array(snp_data[::,req[0]], dtype = float) >= req[1] for req in exclusion_freq], axis = 0, dtype = bool)
                product_mask = product_mask * np.product([np.array(snp_data[::,req[0]], dtype = float) <= req[2] for req in exclusion_freq], axis = 0, dtype = bool)
            else:
                product_mask = True
            exclude += equals_mask * product_mask
    exclude = np.bitwise_not(exclude)
    snp_data = snp_data[exclude]
    snp_categories = np.array(np.zeros(len(snp_data), dtype = int) - 1, dtype = str)
    for category in range(len(classification_scheme)):
        snp_category = np.zeros(len(snp_data), dtype = bool)
        for option in classification_scheme[category]:
            classification_equals = []
            classification_freq = []
            for req in option:
                if req[0] in freq_columns:
                    classification_freq.append(req)
                else:
                    classification_equals.append(req)
            equals_mask = np.product([snp_data[::,req[0]] == str(req[1]) for req in classification_equals], axis = 0, dtype = bool)
            if len(classification_freq) > 0:
                product_mask = np.product([np.array(snp_data[::,req[0]], dtype = float) >= req[1] for req in classification_freq], axis = 0, dtype = bool)
                product_mask = product_mask * np.product([np.array(snp_data[::,req[0]], dtype = float) <= req[2] for req in classification_freq], axis = 0, dtype = bool)
            else:
                product_mask = True
            snp_category += equals_mask * product_mask
        if check_conflicts == True:
            if np.sum(snp_categories[snp_category] == '-1') != np.sum(snp_category):
                raise RuntimeError("%d SNPs of category %d are identified in multiple categories! Check classifier." %(np.abs(np.sum(snp_categories[snp_category] == '-1') - np.sum(snp_category)), category))
        snp_categories[snp_category] = chr(48 + category)
    snp_positions = snp_data[snp_categories != '-1'][::,0]
    snp_categories = np.array(snp_categories[snp_categories != '-1'], dtype = str)
    f_out_type = open(outfile_base + '_types', 'wb')
    f_out_type.write(''.join(snp_categories) + '\n')
    f_out_type.close()
    f_out_pos = open(outfile_base + '_pos', 'wb')
    for pos in snp_positions:
        f_out_pos.write(pos + '\n')
    f_out_pos.close()
    return len(snp_positions), snp_positions[0], snp_positions[-1]
    
