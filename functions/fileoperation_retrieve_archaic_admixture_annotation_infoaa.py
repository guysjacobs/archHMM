###This is a simple program to annotate SNPs on a chromosome with archaic allele sharing information.
###These annotations are then used to build the input of observations used in HMM inference.
###Annotating like this gives flexibility to look for clusters of alleles of different types.

import numpy as np
import gzip
import re

def fileoperation_retrieve_archaic_admixture_annotation_infoaa(infile, outfile, target_name, reference = 'Denisovan', nean_name = '', deni_name = '', pop_list = []):
    """
    Convert each SNP into an annotation based on whether it is shared with Denisovan or Neanderthal, whether it is anc or derived in the target individual, what the maf and daf are in an optionally specified population.

    This classification can then be interpretted as a series of observed 0s and 1s, which can then be classified by the HMM into chunks.
    
    The output order is:
    'POS', 'SNP_A', 'SNP_B', 'SNP_ANC', 'SHARED_A', 'SHARED_B', 'TYPE_A', 'TYPE_B', 'ARCHAIC_HET', 'ARCHAIC_QUALITY', 'ARCHAIC_SHARED', 'ALTARCHAIC_HET',
    where SNP_A, SNP_B, SNP_ANC are 0/1;
    SHARED_A and SHARED_B are 0 if different from target archaic, 1 if the same and -1 if unknown;
    the TYPE_s are -1 (missing) or 0-5:
    0,1,2 = same as archaic target, with ind target being anc, der or unknown
    3,4,5 = different from archaic target, ind target being anc, der or unknown
    ARCHAIC_HET is 0 if hom (most), 1 if het,-1 if missing;
    ARCHAIC_QUALITY is 0 if bad and 1 if good, -1 if missing;
    ARCHAIC_SHARED is 0 if not shared, 1 if shared and -1 if unknown (one or both archaics is 'NN');
    ALTARCHAIC_HET is 0 if hom (most), 1 if het,-1 if missing;
    with optional additional columns POP_MAF and POP_DAF if a population is given.
    """
    open_infile = gzip.open if infile[-3:] == '.gz' else open
    snp_locs = []
    snp_value = [] #Simply the value of the SNP. Useful as we might use whether this is low or high in humans overall as a proxy for ancestral/derived to extend coverage to regions with no Chimp information
    snp_shared = [] #0 if not found in archaic, 1 if found in archaic; -1 if unknown
    snp_type = [] #0 if [ID(test,archaic), test is ancestral], 1 if [ID(test,archaic), test is derived], 2 if [ID(test,archaic), unknown]; 3 if [!ID(test,archaic), test is ancestral], 4 if [!ID(test,archaic), test is derived], 5 if [!ID(test,archaic), unknown]; -1 if unknown
    snp_archaic_het = [] #0 if hom, 1 if het; -1 if missing
    snp_archaic_quality = [] #0 if bad, 1 if good; -1 if missing
    snp_altarchaic_het = [] #0 if hom, 1 if het; -1 if missing
    snp_maf = []
    snp_daf = [] #nan if missing
    snp_archaic_shared = [] #0 if not, 1 if true
    snp_anc = []
    use_pop = False
    snps_read = 0
    with open_infile(infile, 'rb') as f_in:
        for line in f_in:
            if line[0:2] == '##':
                pass
            elif line[0] == '#':
                #headers
                headers = line[0:-1].split('\t')
                idx_target, idx_deni, idx_nean, idx_anc = None, None, None, None
                idx_pop = []
                for i in range(len(headers)):
                    if headers[i] == target_name:
                        idx_target = i if idx_target is None else -1
                    if headers[i] == nean_name:
                        idx_nean = i if idx_nean is None else -1
                    if headers[i] == deni_name:
                        idx_deni = i if idx_deni is None else -1
                    if headers[i] == 'INFO': #ancestral allele is now in the info field
                        idx_anc = i if idx_anc is None else -1
                    if headers[i] in pop_list:
                        idx_pop.append(i)
                if idx_target == -1 or idx_nean == -1 or idx_deni == -1 or idx_anc == -1:
                    print "Multiple individuals with identical names for: %s" %(' '.join([['target', 'neanderthal', 'denisovan', 'ancestral'][i] if [idx_target, idx_nean, idx_deni, idx_anc][i] == -1 else '' for i in range(4)]))
                    raise RuntimeError("Multiple individuals with identical names for: %s" %(' '.join([['target', 'neanderthal', 'denisovan', 'ancestral'][i] if [idx_target, idx_nean, idx_deni, idx_anc][i] == -1 else '' for i in range(4)])))
                if idx_target == None:
                    raise RuntimeError("Could not find target individual with name %s in VCF." %(target_name))
                #print "Retrieved IDx for [%s] and %d population members of %d specified" %(' '.join([['target', 'neanderthal', 'denisovan', 'ancestral'][i] if [idx_target, idx_nean, idx_deni, idx_anc][i] is not None else '' for i in range(4)]), len(idx_pop), len(pop_list))
                comp_idx = idx_nean if reference == 'Neanderthal' else idx_deni if reference == 'Denisovan' else None
                use_pop = True if len(idx_pop) > 0 else False
            else:
                #Retrieve SNP information
                split_line = line[0:-1].split('\t')
                target_snp = split_line[idx_target].split('|')[0:2] #Requires phased. Requires diploid. Expected.
                comp_snp = re.split("[|/]", split_line[comp_idx])[0:2] #Ignores phased. Requires diploid. Expected.
                ref_snp = split_line[3].upper()
                alt_snp = split_line[4].upper()
                anc_snp = split_line[idx_anc].split('AA=')[1][0].upper()
                #Keeping anc_snp consistent with previous notation
                snp_ancestral_inconsistent = False
                if anc_snp == ref_snp:
                    anc_snp = ['0', '0']
                elif anc_snp == alt_snp:
                    anc_snp = ['1', '1']
                elif anc_snp == '.':
                    anc_snp = ['.', '.']
                else:
                    snp_ancestral_inconsistent = True
                snp_locs.append(split_line[1])
                snp_archaic_missing = True if (comp_snp[0] == '.' or comp_snp[1] == '.') else False
                snp_target_missing = True if (target_snp[0] == '.' or target_snp[1] == '.') else False
                if snp_archaic_missing == True or snp_target_missing == True or snp_ancestral_inconsistent == True:
                    if use_pop == True:
                        snp_maf.append('nan')
                        snp_daf.append('nan')
                    snp_shared.append(['nan', 'nan'])
                    snp_type.append(['nan', 'nan'])
                    snp_archaic_het.append('nan')
                    snp_archaic_quality.append('nan')
                    snp_archaic_shared.append('nan')
                    snp_altarchaic_het.append('nan')
                    snp_value.append(['nan', 'nan'])
                    snp_anc.append('nan')
                else:
                    deni_snp = re.split("[|/]", split_line[idx_deni])[0:2] if idx_deni not in [None, -1] else ['N','N'] #Ignores phased. Requires diploid.
                    nean_snp = re.split("[|/]", split_line[idx_nean])[0:2] if idx_nean not in [None, -1] else ['N','N'] #Ignores phased. Requires diploid.
                    altarch_snp = deni_snp if reference == 'Neanderthal' else nean_snp if reference == 'Denisovan' else None
                    snp_value.append([target_snp[0], target_snp[1]])
                    if use_pop == True:
                        num_0 = np.sum([split_line[i].split('|')[0:2].count('0') for i in idx_pop])
                        num_1 = np.sum([split_line[i].split('|')[0:2].count('1') for i in idx_pop])
                        snp_maf.append('%.6f' %(min(num_0, num_1) / float(num_0 + num_1)))
                        snp_daf.append('nan' if anc_snp[0] == '.' else '%.6f' %(num_0 / float(num_0 + num_1)) if anc_snp[0] == '1' else '%.6f' %(num_1 / float(num_0 + num_1)) if anc_snp[0] == '0' else 'nan')
                    snpA_shared = '1' if target_snp[0] in comp_snp else '0'
                    snpB_shared = '1' if target_snp[1] in comp_snp else '0'
                    snp_shared.append([snpA_shared, snpB_shared])
                    if snpA_shared == '1':
                        snpA_type = '2' if anc_snp[0] == '.' else '0' if target_snp[0] == anc_snp[0] else '1'
                    else:
                        snpA_type = '5' if anc_snp[0] == '.' else '3' if target_snp[0] == anc_snp[0] else '4'
                    if snpB_shared == '1':
                        snpB_type = '2' if anc_snp[0] == '.' else '0' if target_snp[1] == anc_snp[0] else '1'
                    else:
                        snpB_type = '5' if anc_snp[0] == '.' else '3' if target_snp[1] == anc_snp[0] else '4'
                    snp_type.append([snpA_type, snpB_type])
                    snp_archaic_het.append('0' if comp_snp[0] == comp_snp[1] else '1')
                    snp_altarchaic_het.append('0' if altarch_snp[0] == altarch_snp[1] else '1')
                    archaic_quality = split_line[7].split('*')
                    if len(archaic_quality) == 1:
                        archaic_quality = '1'
                    else:
                        archaic_quality = archaic_quality[0].split(',')[-1]
                        if reference == 'Denisovan' and archaic_quality[0:3] == 'AD=':
                            archaic_quality = '0'
                        elif reference == 'Neanderthal' and archaic_quality[0:3] == 'AN=':
                            archaic_quality = '0'
                        else:
                            archaic_quality = '1'
                    snp_archaic_quality.append(archaic_quality)
                    snp_archaic_shared.append('-1' if (deni_snp[0] == 'N' and deni_snp[1] == 'N') or (nean_snp[0] == 'N' and nean_snp[1] == 'N')
                                              else '1' if (deni_snp[0] in nean_snp) or (deni_snp[1] in nean_snp)
                                              else '0') #Bugfix V2: new option -1 if either SNP is totally unknown
                    snp_anc.append(anc_snp[0])
                    if snps_read % 50000 == 1:
                        print "pos: ", split_line[1], "\t\tsnp#: ", snps_read
                    snps_read += 1
    print "Writing out file %s with %d snps" %(outfile, snps_read)
    out_open = gzip.open if outfile[-3:] == '.gz' else open
    f_out = out_open(outfile, 'w')
    headers_to_write = ['POS', 'SNP_A', 'SNP_B', 'SNP_ANC', 'SHARED_A', 'SHARED_B', 'TYPE_A', 'TYPE_B', 'ARCHAIC_HET', 'ARCHAIC_QUALITY', 'ARCHAIC_SHARED',  'ALTARCHAIC_HET']
    if use_pop == True:
        headers_to_write.append('POP_MAF')
        headers_to_write.append('POP_DAF')
    f_out.write('\t'.join(headers_to_write) + '\n')
    for snp in range(len(snp_locs)):
        line_to_write = [snp_locs[snp], snp_value[snp][0], snp_value[snp][1], snp_anc[snp], snp_shared[snp][0], snp_shared[snp][1], snp_type[snp][0], snp_type[snp][1], snp_archaic_het[snp], snp_archaic_quality[snp], snp_archaic_shared[snp], snp_altarchaic_het[snp]]
        if use_pop == True:
            line_to_write.append(snp_maf[snp])
            line_to_write.append(snp_daf[snp])
        f_out.write('\t'.join(line_to_write) + '\n')
    f_out.close()
    return None

