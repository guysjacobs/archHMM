# archHMM
This is a HMM to estimate Denisovan/Neanderthal ancestry haplotypes in a target individual, as in the IGDP publication.

The implemetation was initially designed for flexibility (esp. detailed initial annotation and the classification method), allowing the detection of blocks conditioned on multiple archaic hominins for example, or without conditioning on a human outgroup. However, for ease of use this version just includes the HMM described in Jacobs et al 2018.

The sequence of operations is:
a) annotate each SNP based on its allele in target individual, archaic, outgroup etc. (./functions/fileoperation_retrieve_archaic_admixture_annotation_infoaa.py)
b) classify each SNP based its annotation into a string of 1s and 0s (./funcions/fileoperation_classify_archaic_annotation_as_string_array.py)
c) call the Viterbi algorithm on this string, using EM to estimate transition and emission probs (./funcions/calculate_viterbiHMM_binary_archaic_derived_sharing.py)
d) convert the output to a BED file, indicating the span of inferred haplotypes, where they start and end, but not the specific SNPs that contribute.

It runs on one core and is quite slow (e.g. ~one hour per genome). It generates a fair few temporary files during calculations.

The program can be can be called on a VCF or gzipped VCF file, or files. The VCF file must have the ancestral allele in the INFO field as "AA=x", and it must incude at least one archaic genome. A file listing the individuals from the VCF who are the human outgroup populatino, one individual per line, is also needed. If a recombination map is provided then transition probabilities are rates per cM, if not then transition probabilities are probabilities per SNP (and hence vary with genetic diversity, such that even a uniform genetic map is preferable when local recombination rate variation is unknown). Command lines can be adapted from the command given in the /sim/ directory.

The ./sim/ folder contains a simulated 10Mb chromosome following the Malaspinas et al 2016 (https://doi.org/10.1038/nature18299) model of Denisovan introgression into Australians. The command to obtain the output is given in ./sim/\_command.txt . The simulated data was generated using msprime (https://github.com/tskit-dev/msprime, doi: 10.1371/journal.pcbi.1004842) and consist of 35 Africans, 5 Australians and a Denisovan and Neanderthal individual. The actual Denisovan introgressing tracts are recorded in ./sim/sim_35af5au_0.BED.gz , and can be compared to the output. Note that simulated blocks are often mixtures of segments with easy-to-detect coalescent histories - e.g. (Human, (Neanderthal, (Denisovan, Chunk))) - and other coalescent histories, including segments that are unlikely to be identified as introgressed (e.g. ((Human, Chunk), (Neanderthal, Denisovan)). Nevertheless, many blocks are correctly identified, and incorrectly identified blocks are often short or, if longer, caused by Neanderthal introgression.
