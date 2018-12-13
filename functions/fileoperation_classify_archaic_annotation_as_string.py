###This is a file operation that converts summarised SNP data into a string
###of 1s and 0s that can be analysed by the Viterbi algorithm.

###It takes in a file assigning a type to SNPs, discards irrelevant ones, and
###classifies relevant ones according to it's instructions.
###It writes out the SNPs as a single string and their positions as a key file.

###GSJ 21/09/2017

import numpy as np
import copy

def fileoperation_classify_archaic_annotation_as_string(infile, outfile_base, classification_scheme, exclusion_scheme = None, check_conflicts = True):
    """
    Takes in a file assigning a type to SNPs, discards irrelevant ones, and classifies relevant ones according to instructions in classification_scheme.
    classification_scheme is a list of conditions, with [[conditions_0], [conditions_1], ..., [conditions_N]] indicating how each SNP should be classified.
    conditions_0 is a nested list of alternative AND conditions:
    [[[column, value] AND [column, value]] OR [[column, value] AND [column, value]] etc.
    this should be slow but effective.
    The program will check_conflicts by default, but will be faster if this is disabled.
    There is also the option of suggesting an exclusion scheme, which is applied to the snp first and may speed things up by filtering out snps that are not required
    """
    open_f = gzip.open if infile[-3:] == '.gz' else open
    snp_pos = []
    snp_types = []
    conditions_list = []
    #The first step is to convert the classification system to lambda functions. This won't be especially fast...!
    if len(classification_scheme) > 9:
        raise RuntimeError("Only implemented to work with types 0 to 9.")
    #OLD VERSION: classifier_function = lambda snp, classification, conditions: True if np.sum([np.sum([snp[classification[i][j][0]] == str(classification[i][j][1]) for j in range(len(classification[i])))] for i in range(len(classification))], 1) == conditions) > 0 else False
    classifier_function = lambda snp, classification, conditions: True if np.sum([np.sum([snp[classification[i][j][0]] == str(classification[i][j][1]) for j in range(len(classification[i]))]) for i in range(len(classification))] == conditions) > 0 else False
    for classification in classification_scheme:
        conditions_list.append(np.array([len(i) for i in classification]))
        #classifier_list.append(lambda snp, classification, conditions: True if np.sum(np.sum([[snp[classification[i][j][0]] == str(classification[i][j][1]) for j in range(len(classification[i]))] for i in range(len(classification))], 1) == conditions) > 0 else False)
    if exclusion_scheme is None:
        exclusion_conditions = None
        exclusion_lambda = lambda snp, exclusion_scheme, conditions : False
    else:
        exclusion_conditions = np.array([len(i) for i in exclusion_scheme])
        #OLD VERSION exclusion_lambda = lambda snp, exclusion_scheme, exclusion_conditions : True if np.sum(np.sum([[snp[exclusion_scheme[i][j][0]] == str(exclusion_scheme[i][j][1]) for j in range(len(exclusion_scheme[i]))] for i in range(len(exclusion_scheme))], 1) == exclusion_conditions) > 0 else False
        exclusion_lambda = lambda snp, exclusion_scheme, exclusion_conditions: True if np.sum([np.sum([snp[exclusion_scheme[i][j][0]] == str(exclusion_scheme[i][j][1]) for j in range(len(exclusion_scheme[i]))]) for i in range(len(exclusion_scheme))] == exclusion_conditions) > 0 else False
    with open_f(infile, 'rb') as f_in:
        for line in f_in:
            split_line = line[0:-1].split()
            #print split_line
            if exclusion_lambda(split_line, exclusion_scheme, exclusion_conditions) == True:
                pass
            else:
                snp_classified = False
                for classifier in range(len(classification_scheme)):
                    #print snp_classified, classifier_function(split_line, classification_scheme[classifier], conditions_list[classifier])
                    if classifier_function(split_line, classification_scheme[classifier], conditions_list[classifier]) == True:
                        if check_conflicts == False:
                            snp_classified = True
                            snp_type = classifier
                            break
                        else:
                            if snp_classified == True:
                                raise RuntimeError("SNP %s is identified in multiple categories! Check classifier." %(split_line[0]))
                            else:
                                snp_classified = True
                                snp_type = classifier
                    else:
                        pass
                if snp_classified == True:
                    snp_pos.append(split_line[0])
                    snp_types.append(str(snp_type))
    f_out_type = open(outfile_base + '_types', 'wb')
    f_out_type.write(''.join(snp_types) + '\n')
    f_out_type.close()
    f_out_pos = open(outfile_base + '_pos', 'wb')
    for pos in snp_pos:
        f_out_pos.write(pos + '\n')
    f_out_pos.close()
    return None
