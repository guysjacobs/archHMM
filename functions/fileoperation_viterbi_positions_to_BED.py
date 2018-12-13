###This fileoperation combines the output of the Viterbi HMM into a single
###BED file.

###The first few columns of the BED file format are:
#chrom chromStart chromEnd   name  score
#chr7  127471196  127472363  Pos1  0   ...

###I'm going to use chrom, chromStart, chromEnd as usual.
###I'm going to use name to give the type in the HMM (0 or 1 usually) and score
###to indicate how many SNPs support the region.
###I make no statement about what is happening between regions.

###The track is used to indicate the individual and note the meaning of the
###calculation / name.

import gzip
import numpy as np

def fileoperation_viterbi_positions_to_BED(filebase, outfile, type_suffix = 'types_viterbi', chroms = [], track_name = '', track_description = ''):
    #Chroms should be strings
    open_out = gzip.open if outfile[-3:] == '.gz' else open
    with open_out(outfile, 'wb') as f_out:
        f_out.write(' '.join(['track', 'name="%s"' %(track_name), 'description="%s"' %(track_description)]) + '\n')
        for chrom in chroms:
            if type(chrom) is not str:
                raise RuntimeError("Chromosomes in fileoperation_viterbi_positions_to_BED should be given as strings")
            chrom_file_pos = filebase %(chrom) + '_pos'
            chrom_file_types = filebase %(chrom) + type_suffix
            snp_positions = []
            f_chrom_pos = open(chrom_file_pos, 'rb')
            for line in f_chrom_pos:
                snp_positions.append(line[0:-1])
            f_chrom_pos.close()
            f_chrom_types = open(chrom_file_types, 'rb')
            snp_types = np.fromstring(f_chrom_types.readline()[0:-1],dtype='S1') ##LIMITED TO 9 TYPES
            f_chrom_types.close()
            snp_types = np.array(snp_types, dtype = int)
            curr_window_start = None
            curr_window_pos = None
            curr_window_type = None
            num_snps = 0
            for snp in range(len(snp_positions)):
                if curr_window_start is None:
                    curr_window_start = snp_positions[snp]
                    curr_window_pos = snp_positions[snp]
                    curr_window_type = snp_types[snp]
                    num_snps = 1
                elif snp_types[snp] == curr_window_type:
                    #Add and continue
                    curr_window_pos = snp_positions[snp]
                    num_snps += 1
                else:
                    #Write and reset
                    f_out.write('\t'.join([chrom, curr_window_start, '%d' %(int(curr_window_pos)+1), "AN1=%d" %(curr_window_type), "%d" %(num_snps)]) + '\n')
                    curr_window_start = snp_positions[snp]
                    curr_window_pos = snp_positions[snp]
                    curr_window_type = snp_types[snp]
                    num_snps = 1
            f_out.write('\t'.join([chrom, curr_window_start, '%d' %(int(curr_window_pos)+1), "AN1=%d" %(curr_window_type), "%d" %(num_snps)]) + '\n')
    return None
