###This is a program to estimate Denisovan/Neanderthal ancestry haplotypes in a target
###individual.

###fileoperation_retrieve_archaic_admixture_annotation_infoaa saves lots of information about
###each SNP - sharing pattern vs each archaic, MAF and DAF in a population for
###conditioning, sharing pattern between archaics, anc/der etc.

###For now, I:
###a) archaic-annotate sites
###b) classify sites to a string of 1s and 0s
###c) call the viterbi program on this string
###d) convert the output to a BED file, indicating the span of inferred haplotypes,
###where they start and end, but not the specific SNPs that contribute.

###Things to potentially add to increase functionality:
###Forward-Backward agorithm
###More complex models - more categories or types of transition
###Different weight to different types of allele sharing, ideally based on different population trees.

###There are initial parameters and a penalty function for the Viterbi. Higher penalty means fewer
###tracts, which will either mean that separated archaic SNPs are more likely to be linked up
###or will lead to a relaxed criteria for how many archaics in a 0-tract such that they are no
###longer linked. Probably the latter.

###NOTE: This is a relatively slow program that works better for small groups of individuals.
###If there are 400 individuals, it needs to read in and analyse each chromosome for each
###individual. A full genome for a single individual can easily take an hour.
###I intend to write an alternative program that is much faster.

###GSJ written September 2017 (03/10/17) v1.0
###V1.1 GSJ 08/11/17:
##Modified the exclusion schemes to exclude SNPs with missing data in either archiac in binary_archaic_derived_sharing, binary_archaic_derived_sharing_no_derived_divergence and binary_archaic_derived_sharing_no_derived_divergence_joint_divergence
##Added in a functionality whereby chosen parameters do not update during the Viterbi learning iterations. This also involved changes to calculate_viterbiHMM_binary_archaic_derived_sharing.
##Added in a functionality whereby the transition probability is estimated for a set distance, based on the chromosome length and number of SNPs. This involved asking fileoperation_classify_archaic_annotation_as_string_array to return the number of SNPs and start and end of the chromosome.
##Added in a functionality whereby the HMM can be run multiple times, using either multiple specified initial conditions or randomly chosen conditions. This also involved changes to calculate_viterbiHMM_binary_archaic_derived_sharing.
##Rectified the saving of the a single SNP chunk to [start,end] = [x,x+1] in the BED rather than [x,x]. This involved changes to fileoperation_viterbi_positions_to_BED.

###V1.2 GSJ 09/11/17:
##Added in support for fitting the transition rate, and asymmetric rates, in the binary_archaic_derived_sharing HMM. This also involved changes to calculate_viterbiHMM_binary_archaic_derived_sharing .
##Removed the option of estimating the transition probability based on SNP spacing, which lead to very different chunk sizes between chromosomes.
##Added in support for variable recombination rate (in binary_allele_sharing model)
##Adding in support for the six category model (in progress)
##Adding in support for assessing power (in progress)
#The point of this is to keep access to the different functions under one 'roof'
#I am interested in the false and true positive rates in a non-admixed and admixed population
#To do this I need to:
#Generare samples of two archaics, an admixed chromosome and a non-admixed chromosome, using scrm.
#Record the local geneologies. If there is coalescence before the human/archiac MRCA (assuming 1-way
#migration) then this is introgressed. Otherwise, it isn't. Note that actually things are more complex.
#Assign chunks according to the trees.
#Repeat *replicate* number of times.
#Run the HMM according to parameters on each case.
#Record summary statistics about the output - the proportion of admixture that is correctly or
#incorrectly inferred, the average error in the end positions when there is an overlap...

###V1.3 GSJ 22/11/17:
##Added in support for a profiling mode, which just counts out types per chromosome without doing viterbi.
##Added in support for searching for local admixture

###V1.4 GSJ 02/02/2018
##Modified so that it uses the new VCF files, in which the ancestral allele is in the INFO field as AA=.

#TO DO:
#Add in support for the joint estimation of Nean and Deni.
#Add in support for a Rampasassa-centric method.
#The model is that some Rampasassa will have more hobbit. The hobbit is a non-human ancestor.
#I'm looking for SNPs that are higher freq in Ram, 


import sys
import numpy as np
import argparse
import os
import subprocess
import time
import gzip

#sys.path.insert(0, '/media/sf_Dropbox/PythonScripts/') #If this location doesn't exist it shoudldn't matter; but the dependencies need to be in the same directory as this file.
sys.path.insert(0, 'I:/Dropbox/PythonScripts/')
sys.path.insert(0, '/mmg/jacobs/PythonPrograms/functions/')
print sys.path
from fileoperation_classify_archaic_annotation_as_string_array import *
from fileoperation_retrieve_archaic_admixture_annotation_infoaa import *
from fileoperation_log_append_from_list import *
from fileoperation_interpolate_gmap_1KG import *
from calculate_viterbiHMM_binary_archaic_derived_sharing import calculate_viterbiHMM_binary_archaic_derived_sharing, support_calculate_viterbi_probabililty
from fileoperation_viterbi_positions_to_BED import *
from readin_populationFileEOL import *

parser = argparse.ArgumentParser(description='Calculate tracts of archaic introgression based on a Viterbi algorithm and allele sharing, based on a VCF file including both a target individual and archaics.')

parser.add_argument('--infile', metavar='infile', type=str, nargs='*', default = ['None'],
                    help='system location of the VCF input file; multiple e.g. chromosome files can be passed, but if so their respective chromosomes must still be specified using --chromosomes.')
parser.add_argument('--outfile', metavar='outfile', type=str, nargs='?', default = 'None',
                    help='system location of the output file root. This program generates two outputs: _tracts (indicating when different tracts begin and end) and _probs (indicating the fitted archaic frequencies)')

parser.add_argument('--infile_recom', metavar='infile_recom', type=str, nargs='*', default = [],
                    help='system location of a genetic map input file, if available. Multiple e.g. chromosome files can be passed, but if so their respective chromosomes must still be specified using --chromosomes.')

parser.add_argument('--target_individuals', '-targi', dest='target_individuals', type=str, action = 'store', nargs='*', default = [],
                    help='which individuals to calculate archaic admixture tracts for; this is the name in the VCF.')
parser.add_argument('--target_archaic', '-targa', dest='target_archaic', type=str, action = 'store', nargs='?', default = '',
                    help='calculate on the Denisovan or Neanderthal? NB I hope to increase flexibility later, especially re: testing sharing tracts vs humans.')
parser.add_argument('--denisovan_name', '-deni', dest='denisovan_name', type=str, action = 'store', nargs='?', default = '',
                    help='which individual is the Denisovan in the VCF.')
parser.add_argument('--neanderthal_name', '-nean', dest='neanderthal_name', type=str, action = 'store', nargs='?', default = '',
                    help='which individual is the Neanderthal in the VCF.')
#parser.add_argument('--ancestral_name', '-anc', dest='ancestral_name', type=str, action = 'store', nargs='?', default = '',
#                    help='which individual is the Ancestral in the VCF.')
parser.add_argument('--chromosomes', '-c', dest='chromosomes', action = 'store', nargs = '*', default = [],
                    help="list (or range, '1 ... 22') of chromosomes to calculate archaic tracts on. These correspond in order to the files given as infile, or are substituted as integers into a single infile at the '%%d' position.")
parser.add_argument('--log', '-l', dest='log', action = 'store_true', default = False,
                    help="should log files be written? Not implemented yet.")
parser.add_argument('--popfile', '-popfile', dest='popfile', type=str, action = 'store', nargs = '?', default = '',
                    help="population file (one individual name per line) to be used when allele classes are conditioned on population frequency")

parser.add_argument('--viterbi_method', '-vm', dest='viterbi_method', action = 'store', nargs='?', default = 'binary_archaic_derived_sharing',
                    help='the HMM method to be applied. This impacts how the archaic annotation is summarised and the HMM algorithm used. Methods based on "binary_archaic_derived_sharing(...)" use the same algorithm, but classify different types of SNPs as evidence for or against admixture. I hope to implement more complex HMMs in the future, or add subtlety to the annotation method.')
parser.add_argument('--viterbi_transition', '-vt', dest='viterbi_transition', type=float, action = 'store', nargs='*', default = [0.9,0.1,0.1,0.9],
                    help='the transition probability matrix for switching states. Low tranistion probabilities lead to longer haplotypes. The full matrix should be specified i.e. for two hidden states: a, b, c, d is [[a,b],[c,d]] is [[fromAtox], [fromBtox]], and a and d are the probabilities of no transition, and and a should be 1.0 - b and d is 1.0 - c. I now offer an option to fit this parameter using viterbi_fit_transition_rates .')
parser.add_argument('--viterbi_init_parameters', '-vparam', dest='viterbi_init_parameters', type=float, action = 'store', nargs='*', default = [],
                    help='the initial parameters for the viterbi algorithm. The meaning depends on the algorithm chosen. For binary_archaic_derived_sharing, we expect three parameters, the probability of the first state (left most in the genome) being each hidden state (human/archaic) which should not matter, the rate of archaic signal in low-archaic regions and the rate of archaic signal in high-archaic regions.')
parser.add_argument('--viterbi_init_resamples', '-vnuminit', dest = 'viterbi_init_resamples', type = int, action = 'store', nargs = '?', default = '1',
                    help="how many different (random) initial conditions are attempted. If this is > 1 then the HMM will make multiple attempts at optimisation, using randomly sampled initial conditions from uniform distributions specified in viterbi_init_parameters.")
parser.add_argument('--viterbi_init_alternatives', '-valtinit', dest = 'viterbi_init_alternatives', type = int, action = 'store', nargs = '?', default = '1',
                    help="are multiple initial conditions provided? If so, these are given as A1 A2 A3 B1 B2 B3 ... CN CN CN where A,B,C are model inputs and provided in viterbi_init_parameters. viterbi_init_alternatives should be set to N, and indicates how many different conditions are attempted.")
parser.add_argument('--viterbi_fixed_parameters', '-vparamfix', dest='viterbi_fixed_parameters', type=int, action = 'store', nargs='*', default = [],
                    help='the IDxs of the parameters to keep fixed when optimising the chunks. This can be relevant for between-population comparisons or if a single introgression event is assumed. Be aware that this radically changes what the program is telling us. In the case of binary_archaic_derived_sharing and similar, the parameter IDxs are 0 for the initial amount of introgressive SNPs in the low class and 1 for the initial amount of introgressive SNPs in the high class.')
parser.add_argument('--viterbi_fit_transition_rates', '-vfittrans', dest='viterbi_fit_transition_rates', action = 'store_true', default = False,
                    help="fir the transition rates as well as the probabilities. This involves counting transitions along the chromosome after each Viterbi iteration, and updating. Best used when simply searching for the best model for a specific chromosome. May have high sensitivity to initial conditions.")


parser.add_argument('--cleanup_level', '-cl', dest = 'cleanup_level', type = int, action = 'store', nargs = '?', default = '0',
                    help="which files to delete during the calculation. 0 means no files are deleted; 1 deletes the archaicAnnotation file.")
parser.add_argument('--chrom_to_bed', '-chrombed', dest='chrom_to_bed', action = 'store_true', default = True,
                    help="convert the Viterbi output files to a BED file for comparing individuals/ease of manipulation.")

parser.add_argument('--profilemode', '-prof', dest='profilemode', action = 'store_true', default = False,
                    help="run the program in profile mode? In this case, no simulations or HMM is run. Instead, the SNPs are classified according to the viterbi_method annotation, and then the categories are simply counted. This is saved as a per-chromosome and total file.")

parser.add_argument('--powermode', '-pow', dest='powermode', action = 'store_true', default = False,
                    help="run the program in power mode? In this case, simulations are used rather than the input file, and the known output of the simulations is compared to the input.")
parser.add_argument('--sim_program_location', '-simprogloc', dest='sim_program_location', action = 'store', nargs='?', default = 'I:/PhylogenPrograms/scrm/scrm.exe',
                    help="the system location of the simulation program to be used. I use scrm which takes ms style input and returns local trees using the -T tag.")

parser.add_argument('--sim_command_line_override', '-simcommandover', dest='sim_command_line_override', action = 'store', nargs='?', default = '',
                    help="this can be used to apply a custom simulation command line to scrm. Otherwise, I build the command myself based on a simple introgression model.")
parser.add_argument('--sim_replicates', '-simreps', dest='sim_replicates', type = int, action = 'store', nargs = '?', default = '1',
                    help="how many replicates?")
parser.add_argument('--sim_locus_length', '-simlength', dest='sim_locus_length', type = int, action = 'store', nargs = '?', default = '1',
                    help="how many bp is the locus being simulated?")

parser.add_argument('--sim_model', '-simmodel', dest='sim_model', action = 'store', nargs='?', default = 'kuhlwilm2016',
                    help="choose which model to use. The options are kuhlwilm2016 and rogers2017. rogers2017 involves an estimated Denisovan population size equal to the Neanderthal one. Note that I add in a chimp in the modelling, with an (immaterial) population size set to 25000 in approximate corresponance to Prado-Martinez et al 2013.")

"""
parser.add_argument('--sim_mutation_rate', '-simmutrate', dest='sim_mutation_rate', type = int, action = 'store', nargs = '?', default = '1',
                    help="how many mutations per year are assumed? This is used to scale the split times and calculate thetas etc.")
parser.add_argument('--sim_recombination_rate', '-simrecrate', dest='sim_recombination_rate', type = int, action = 'store', nargs = '?', default = '1',
                    help="how many recombinations per site per year are assumed?")

parser.add_argument('--sim_introgression_proportion', '-simintroprop', dest='sim_introgression_proportion', type = float, action = 'store', nargs = '?', default = '0.02',
                    help="how much introgression?")
parser.add_argument('--sim_introgression_time', '-simintrotime', dest='sim_introgression_time', type = float, action = 'store', nargs = '?', default = '0.02',
                    help="when did introgresion happen, in generations ago?")
parser.add_argument('--sim_split_deninean', '-simsplitDN', dest='sim_split_deninean', type = int, action = 'store', nargs = '?', default = '0.02',
                    help="when did Denisovan and Neanderthal split, in generations?")
parser.add_argument('--sim_split_denineanhum', '-simsplitDNH', dest='sim_split_denineanhum', type = int, action = 'store', nargs = '?', default = '0.02',
                    help="when did Human split with Denisovan and Neanderthal, in generations?")
parser.add_argument('--sim_split_denineanhumchimp', '-simsplitDNHC', dest='sim_split_denineanhumchimp', type = int, action = 'store', nargs = '?', default = '0.02',
                    help="when did Chimp split with Human, Denisovan and Neanderthal, in generations? NB that the chimp is needed to guess the ancestral state.")
parser.add_argument('--sim_popsize_hum', '-simpopH', dest='sim_popsize_hum', type = float, action = 'store', nargs = '?', default = '61231',
                    help="haploid population size of Human?")
parser.add_argument('--sim_popsize_deni', '-simpopD', dest='sim_popsize_deni', type = float, action = 'store', nargs = '?', default = '31166',
                    help="haploid population size of Denisovan? This wasn't estimated by Rogers et al, so I may need to use ")
parser.add_argument('--sim_popsize_nean', '-simpopN', dest='sim_popsize_nean', type = float, action = 'store', nargs = '?', default = '31166',
                    help="haploid population size of Neanderthal?")
parser.add_argument('--sim_popsize_chimp', '-simpopC', dest='sim_popsize_chimp', type = float, action = 'store', nargs = '?', default = '10000',
                    help="haploid population size of Chimp? This isn't important as no introgression etc.")
parser.add_argument('--sim_popsize_deninean', '-simpopDN', dest='sim_popsize_deninean', type = float, action = 'store', nargs = '?', default = '637',
                    help="haploid population size of Denisovan/Neanderthal common ancestral population?")
parser.add_argument('--sim_popsize_denineanhum', '-simpopDNH', dest='sim_popsize_denineanhum', type = float, action = 'store', nargs = '?', default = '19195',
                    help="haploid population size of Denisovan/Neanderthal/Human common ancestral population?")
"""


args = parser.parse_args()

##The methods are:
#binary_archaic_derived_sharing
#binary_archaic_derived_sharing_no_derived_divergence
#binary_archaic_derived_sharing_no_derived_divergence_joint_divergence
#six_category_simple_archaic_demography
#fourteen_category_explore_archaic_pop_sharing
"""
##For running from Python shell
args.infile = ["I:/Genetic_Data/Indonesia_Diversity175/multicall3_mask1/phase/VCFsWithArchaic/chr%d_cr99_merged_PIR.phased.v2.haps.AN.AD.AA_as_individuals.479.vcf.gz"]
args.viterbi_method = 'binary_archaic_derived_sharing_racimo2016' #'six_category_simple_archaic_demography' #'binary_archaic_derived_sharing_no_derived_divergence'
short_method = 'bads' if args.viterbi_method == 'binary_archaic_derived_sharing' else 'badsndd' if args.viterbi_method == 'binary_archaic_derived_sharing_no_derived_divergence' else 'badsnddjd' if args.viterbi_method == 'binary_archaic_derived_sharing_no_derived_divergence_joint_divergence' else 'scsad' if args.viterbi_method == 'six_category_simple_archaic_demography' else 'other'
args.outfile = "I:/Genetic_Data/Indonesia_Diversity175/multicall3_mask1/phase/VCFsWithArchaic/Allele_sharing/racimo_test_deni/14catAFracimo_%s_%s_fix060_denisovan" %(short_method, '%s')

args.infile_recom = ["I:/Genetic_Data/1KGenomes/1000GP_Phase3/genetic_map_chr%d_combined_b37.txt"]
args.target_individuals = ['MPI-025', 'LP6005441-DNA_A03']#['SS6004468', 'NIAS21', 'MPI-025']#, 'LP6005441-DNA_E07'] #MPI-025 NIAS21 LP6005441-DNA_E07
#args.target_individuals = ['LP6005441-DNA_D05', 'LP6005441-DNA_C05', 'SS6004469', 'LP6005442-DNA_E10', 'LP6005442-DNA_F10', 'LP6005442-DNA_G09', 'LP6005442-DNA_H09', 'LP6005442-DNA_F01', 'LP6005443-DNA_D02', 'LP6005441-DNA_G03', 'LP6005441-DNA_H03'] #For testing Nean
args.target_individuals = ['MPI-025', 'MPI-030', 'MPI-065', 'MPI-070', 'UV002', 'UV003', 'UV355B', 'UV368B', 'LP6005441-DNA_A03', 'LP6005441-DNA_B03', 'LP6005442-DNA_E10', 'LP6005442-DNA_F10', 'LP6005442-DNA_G03', 'LP6005442-DNA_H03'] #For testing Deni
args.popfile = "I:/Genetic_Data/Indonesia_Diversity175/population_data/sample_lists/list_subset_africa_SSonly.txt" #"I:/Genetic_Data/Indonesia_Diversity175/population_data/sample_lists/list_subset_mainlandseataiwanigorotpapua.txt"
##Uncomment the following to calcualte for all individuals; doing East Asians, Southeast Asians, Island Southeast Asians and Papuans first.
#args.target_individuals = ['ALR03','ALR04','ALR06','ALR11','ALR21','BLI615','TL003','TL004','TL032','TL038','TL097','TL113','TL119','LP6005519-DNA_E06','LP6005519-DNA_F06','LP6005441-DNA_D04','LP6005443-DNA_B01','LP6005592-DNA_D03','SS6004467','LP6005441-DNA_F04','LP6005441-DNA_C05','LP6005441-DNA_D05','SS6004469','LP6005441-DNA_G05','LP6005441-DNA_H05','LP6005441-DNA_C06','LP6005441-DNA_D06','LP6005592-DNA_C02','LP6005443-DNA_C06','LP6005443-DNA_D06','LP6005441-DNA_B07','LP6005443-DNA_E01','LP6005441-DNA_C08','LP6005441-DNA_D08','LP6005441-DNA_B09','LP6005443-DNA_E09','LP6005441-DNA_E09','LP6005441-DNA_F09','LP6005443-DNA_F01','LP6005443-DNA_G01','LP6005441-DNA_D12','LP6005443-DNA_H01','LP6005441-DNA_F12','LP6005443-DNA_A02','LP6005442-DNA_B01','LP6005443-DNA_B02','LP6005442-DNA_D01','LP6005443-DNA_C02','LP6005442-DNA_G01','LP6005442-DNA_H01','BNA01','BNA03-F','BNA05','BNA12-F','BNA14-F','BNA21-F','BNA22-F','BNA25-F','BNA26-F','BNA29-F','BNA40-F','BRE05','BRE06','BRE10','CBL002','CBL006','CBL010','CBL018','CBL019','CBL022','CBL025','CBL027','CBL054','CBL33','CBL34','CBL49','CBL55','RAM003','RAM005','RAM008','RAM022','RAM024','RAM025','RAM027-F','RAM034-F','RAM035-F','RAM036-F','RAM038','RAM039-F','RAM041-F','RAM043-F','RAM045-F','RAM067','RAM087','RAM105','DNG02','DNG05','DNG06','DNG09','DNG17','DNG34','DTT-KEI-006','DTT-KEI-031','FAN-KEI-011','FAN-KEI-024','WAR-KEI-005','WAR-KEI-030','FLT018','FLT020','FLT022','FLT042','FLT060','FLT077','FLT078','MTW007','MTW020','MTW030','MTW057','MTW058','MTW062','MTW071','MTW078','NIAS15','NIAS16','NIAS17','NIAS21','NIAS26','NIAS33','NIAS38','NIAS01','NIAS02','NIAS03','NIAS04','NIAS08','NIAS10','NIAS12','NIAS13','LP6005592-DNA_B02','LP6005592-DNA_H03','SS6004477','SS6004478','UV002','UV003','UV008','UV030','UV031','UV034','UV275B','UV287B','UV291','UV293','UV321B','UV322B','UV344B','UV350','UV355B','UV368B','LP6005441-DNA_A03','LP6005441-DNA_B03','MPI-025','MPI-030','MPI-065','MPI-070','MPI-074','MPI-296','MPI-376','NG02-F','NG06','NG09-F','NG63-F','NG65-F','NG66','NG88-F','PNG001','PNG003','PNG005','PNG031','PNG057','PNG059','PNG064','PNG076','LP6005441-DNA_A10','LP6005441-DNA_B10','LP6005443-DNA_A08','LP6005443-DNA_B08','LP6005443-DNA_C07','LP6005443-DNA_C08','LP6005443-DNA_D07','LP6005443-DNA_D08','LP6005443-DNA_E07','LP6005443-DNA_E08','LP6005443-DNA_F07','LP6005443-DNA_F08','LP6005443-DNA_G07','LP6005443-DNA_H07','SS6004472','EGAN00001279031','EGAN00001279032','EGAN00001279033','EGAN00001279034','EGAN00001279035','EGAN00001279036','EGAN00001279037','EGAN00001279038','EGAN00001279039','EGAN00001279040','EGAN00001279041','EGAN00001279042','EGAN00001279043','EGAN00001279044','EGAN00001279045','EGAN00001279046','EGAN00001279047','EGAN00001279048','EGAN00001279049','EGAN00001279050','EGAN00001279051','EGAN00001279052','EGAN00001279053','EGAN00001279054','EGAN00001279055','LP6005519-DNA_C06','LP6005519-DNA_D06','LP6005519-DNA_A06','LP6005519-DNA_B06','LP6005441-DNA_G03','LP6005441-DNA_H03','LP6005443-DNA_A07','LP6005443-DNA_B07','LP6005442-DNA_C11','LP6005442-DNA_D11','KJG02','KJG04','KJG08','KJG18','KJG25','KJG32','SLW13','SLW36','SLW40','SLW42','SLW51','SLW54','SMT01','SMT03','SMT05','SMT12','SMT22','SMT23','SMT39','LP6005442-DNA_C07','LP6005442-DNA_E07','LP6005443-DNA_G05','ALK-TBR-025','FDT-TBR-004','MKT-TBR-007','SLD-TBR-001','SLD-TBR-022','TMB-TBR-027','LP6005441-DNA_A08','LP6005441-DNA_A11','LP6005441-DNA_B02','LP6005441-DNA_B08','LP6005441-DNA_E07','LP6005441-DNA_F07','LP6005441-DNA_G02','LP6005441-DNA_G08','LP6005441-DNA_H02','LP6005441-DNA_H08','LP6005442-DNA_A02','LP6005442-DNA_A10','LP6005442-DNA_B02','LP6005442-DNA_B10','LP6005442-DNA_D09','LP6005442-DNA_E11','LP6005442-DNA_F09','LP6005442-DNA_F11','LP6005442-DNA_G10','LP6005442-DNA_G11','LP6005442-DNA_H10','LP6005442-DNA_H11','LP6005443-DNA_A01','LP6005443-DNA_B09','LP6005443-DNA_E02','LP6005443-DNA_E06','LP6005443-DNA_F02','LP6005443-DNA_F06','LP6005443-DNA_G02','LP6005443-DNA_G08','LP6005443-DNA_H08','LP6005592-DNA_C03','LP6005592-DNA_C05','LP6005619-DNA_B01','LP6005619-DNA_C01','LP6005677-DNA_D03','LP6005677-DNA_G01','SS6004470','SS6004471','SS6004473','SS6004475','SS6004480','LP6005441-DNA_A12','LP6005441-DNA_B04','LP6005441-DNA_B12','LP6005441-DNA_E10','LP6005441-DNA_F10','LP6005441-DNA_G06','LP6005441-DNA_G07','LP6005441-DNA_H06','LP6005441-DNA_H07','LP6005442-DNA_A07','LP6005442-DNA_B07','LP6005442-DNA_H06','LP6005443-DNA_A12','LP6005443-DNA_E11','LP6005443-DNA_F05','LP6005443-DNA_F11','LP6005443-DNA_G11','LP6005443-DNA_H11','LP6005519-DNA_A03','LP6005519-DNA_B03','LP6005519-DNA_D01','LP6005519-DNA_G02','LP6005677-DNA_D01','LP6005677-DNA_E01','LP6005677-DNA_F01','SS6004476','SS6004479','LP6005441-DNA_E08','LP6005441-DNA_F08','LP6005442-DNA_E12','LP6005442-DNA_F01','LP6005442-DNA_F02','LP6005442-DNA_F12','LP6005442-DNA_G12','LP6005442-DNA_H12','LP6005443-DNA_A03','LP6005443-DNA_B03','LP6005443-DNA_B04','LP6005443-DNA_C03','LP6005443-DNA_C04','LP6005443-DNA_D02','LP6005443-DNA_D03','LP6005443-DNA_D04','LP6005443-DNA_D05','LP6005443-DNA_E05','LP6005443-DNA_F03','LP6005443-DNA_F04','LP6005443-DNA_G03','LP6005443-DNA_G04','LP6005443-DNA_H02','LP6005443-DNA_H03','LP6005592-DNA_F03','LP6005677-DNA_A02','LP6005677-DNA_B02','LP6005441-DNA_C01','LP6005441-DNA_C03','LP6005441-DNA_C07','LP6005441-DNA_C10','LP6005441-DNA_D01','LP6005441-DNA_D03','LP6005441-DNA_D07','LP6005441-DNA_D10','LP6005441-DNA_E03','LP6005441-DNA_E05','LP6005441-DNA_E06','LP6005441-DNA_F03','LP6005441-DNA_F05','LP6005441-DNA_F06','LP6005441-DNA_G11','LP6005441-DNA_H11','LP6005442-DNA_A12','LP6005442-DNA_B12','LP6005442-DNA_G09','LP6005442-DNA_H09','LP6005443-DNA_C09','LP6005443-DNA_D09','LP6005519-DNA_A04','LP6005519-DNA_A05','LP6005519-DNA_B04','LP6005519-DNA_B05','LP6005519-DNA_C04','LP6005519-DNA_C05','LP6005519-DNA_D04','LP6005519-DNA_D05','LP6005519-DNA_E04','LP6005519-DNA_E05','LP6005519-DNA_F04','LP6005519-DNA_G03','LP6005519-DNA_G04','LP6005519-DNA_H03','LP6005519-DNA_H04','LP6005592-DNA_A04','LP6005592-DNA_B04','LP6005592-DNA_E05','LP6005592-DNA_F05','LP6005592-DNA_G05','LP6005619-DNA_A01','LP6007071-DNA_A01','LP6007072-DNA_A01','LP6007073-DNA_A01','LP6007074-DNA_A01','LP6007075-DNA_A01','LP6007076-DNA_A01','LP6005441-DNA_A01','LP6005441-DNA_A05','LP6005441-DNA_A06','LP6005441-DNA_B01','LP6005441-DNA_B05','LP6005441-DNA_B06','LP6005441-DNA_C02','LP6005441-DNA_C09','LP6005441-DNA_C11','LP6005441-DNA_D02','LP6005441-DNA_D09','LP6005441-DNA_D11','LP6005441-DNA_E02','LP6005441-DNA_F02','LP6005441-DNA_G04','LP6005441-DNA_G09','LP6005441-DNA_G10','LP6005441-DNA_G12','LP6005441-DNA_H09','LP6005441-DNA_H10','LP6005441-DNA_H12','LP6005442-DNA_A03','LP6005442-DNA_A04','LP6005442-DNA_A08','LP6005442-DNA_A11','LP6005442-DNA_B03','LP6005442-DNA_B04','LP6005442-DNA_B08','LP6005442-DNA_B11','LP6005442-DNA_C02','LP6005442-DNA_C04','LP6005442-DNA_C10','LP6005442-DNA_D02','LP6005442-DNA_D03','LP6005442-DNA_D08','LP6005442-DNA_D10','LP6005442-DNA_E04','LP6005442-DNA_E10','LP6005442-DNA_F04','LP6005442-DNA_F10','LP6005442-DNA_G02','LP6005442-DNA_G03','LP6005442-DNA_G04','LP6005442-DNA_G07','LP6005442-DNA_H03','LP6005442-DNA_H04','LP6005443-DNA_A06','LP6005443-DNA_B06','LP6005443-DNA_B10','LP6005443-DNA_D01','LP6005443-DNA_E10','LP6005443-DNA_F10','LP6005443-DNA_H05','LP6005519-DNA_C03','LP6005519-DNA_D03','LP6005519-DNA_F03','LP6005592-DNA_A02','LP6005592-DNA_B01','LP6005592-DNA_B03','LP6005592-DNA_C01','LP6005592-DNA_D01','LP6005592-DNA_D04','LP6005592-DNA_E01','LP6005592-DNA_E02','LP6005592-DNA_F01','LP6005592-DNA_G01','LP6005592-DNA_G03','LP6005592-DNA_H01','LP6005677-DNA_A03','LP6005677-DNA_B01','LP6005677-DNA_C03','LP6007068-DNA_A01','LP6007069-DNA_A01','SS6004468','SS6004474']
args.target_archaic = 'Denisovan' #'Neanderthal'
args.denisovan_name = 'DENI'
args.neanderthal_name = 'NEAN'
args.chromosomes = ['1','...','22']
#args.viterbi_transition = [0.9999, 0.0001, 0.1, 0.9] #0.1 implies an average haplotype is 10 SNPs, which might be e.g. 1600*10 = 16kb. Others have used a switch every 0.00005M = 0.005cM = 5000bp...
args.viterbi_transition = [100.0, 100.0, 100.0, 100.0] #0.1 implies an average haplotype is 10 SNPs, which might be e.g. 1600*10 = 16kb. Others have used a switch every 0.00005M = 0.005cM = 5000bp...
args.viterbi_init_parameters = [0.01, 0.1, 0.60] #For binary_archaic_derived_sharing this is the probability of the first (left-most) SNP being in each archaic introgression category (unimportant); the rate of archaic signal in low-archaic regions; and the rate of archaic signal in high-archaic regions.
args.viterbi_fit_transition_rates = True

#args.viterbi_init_alternatives = 1
#args.viterbi_init_resamples = 1
args.viterbi_fixed_parameters = [1]

#args.sim_model = 'malaspinas2016'
#args.powermode = False
#args.sim_locus_length = 1000000

args.profilemode = False

args.log = True
args.cleanup_level = 2 #This should be 0 when testing; 1 deletes the actual SNP summary information while 2 deletes the Viterbi strings as well (only use when saving as BED!)
args.chrom_to_bed = True
"""
print 'Input file(s): ',  ' '.join(['%s' %(infile) for infile in args.infile])
print 'Output file: ', args.outfile
print 'Target individual list: ', ' '.join(['%s' %(ind) for ind in args.target_individuals])
#print 'Ancestral population list: ', ' '.join(['%s' %(pop) for pop in args.ancestral_population_list])
print 'Chromosomes used: ', ' '.join(['%s' %(chrom) for chrom in args.chromosomes])
print "Initial conditions: ", args.viterbi_init_parameters
print "Transition inital conditions: ", args.viterbi_transition

if (args.infile[0] == ["None"] and args.powermode is False) or args.outfile == "None":
    print "Please specify input and output files."

elif (len(args.target_individuals) == 0 and args.powermode is False) or (args.target_archaic != 'Neanderthal' and args.target_archaic != 'Denisovan'):
    print "Please specify at least one target individual and either Denisovan or Neanderthal as the target individual (%s)." %(args.target_archaic)

elif args.neanderthal_name == '' or args.denisovan_name == '':
    print "Please specify the name of the Neanderthal and Denisovan individuals."

elif len(args.target_individuals) > 1 and args.outfile.count('%s') != 1:
    print "When multiple individuals are specified, the outfile must contain a single substitution position '%s'"

elif len(args.infile) > 1 and not (len(args.infile) == len(args.chromosomes) or (len(args.chromosomes) == 3 and args.chromosomes[1] == '...' and args.chromosomes[2] - args.chromosomes[0] == len(args.infile) - 1)):
    print "Multiple VCF files contributing to the calculation need to have their chromosomes specified in order using --chromosomes."

elif len(args.infile_recom) > 0 and (len(args.infile_recom) != len(args.infile)):
    print "If a genetic map is used, provide the same number of maps as VCFs, or, if a '%d' is used in one, use in the other also"

elif len(args.infile) == 1 and len(args.chromosomes) == 0:
    print "When a single infile is specified, at least one chromosome should be specified using --chromosomes.  Setting this to a dummy variable (e.g. 0, -1, 99) should work but isn't tested; be careful if program is extended to e.g. cut to exons."

elif len(args.chromosomes) > 1 and args.infile[0].count('%d') != 1:
    print "When a single infile is specified, and multiple chromosomes are specified using --chromosomes the infile must contain a single substitution position '%d'"

elif args.viterbi_init_alternatives > 1 and args.viterbi_init_resamples > 1:
    print "When using multiple initial conditions, you can only use random sampled conditions (--viterbi_init_resamples > 1) or multiple input conditions (--viterbi_init_alternatives > 1), not both"

elif (args.infile_recom) == 0 and np.sum(args.viterbi_transition) != np.sqrt(len(args.viterbi_transition)):
    print "viterbi_transition is the transition matrix, which are probabilities if no recombination map is supplied and so should sum to K (number of hidden states)."

elif args.profilemode is True and args.powermode is True:
    print "Only one of profilemode and powermode can be True. Use profilemode to simply count the number of SNPs in different categories (following viterbi_method) and powermode to assess the power of the HMM via simulations"

elif args.viterbi_method is 'fourteen_category_explore_archaic_pop_sharing' and args.popfile is '':
    print "The method fourteen_category_explore_archaic_pop_sharing incorporates an out or in-group population, and assesses sharing between the target individual, Nean, Deni and that population. Provide a population using --popfreq "

else:
    #If the program is doing a power calculation, generate the command line and run the simulation.
    if args.powermode is True:
        #In powermode!
        if args.sim_command_line_override != '':
            simulation_command = args.sim_command_line_override
        elif args.sim_model == 'rogers2017':
            pass
            """Ne_N, Ne_H, Ne_DN, Ne_DNH = [61231, 31166, 637, 19195]
            T_N, T_DN, T_DNH = [1897, 25660, 25920] #generations, divergence times
            S_D, S_N = [1734.5, 1760] #generations, sampling times
            Ne_D = 31166 #THIS IS MADE UP, set as the same as Neanderthal
            theta_0 = 4.0 * args.sim_popsize_Hum * args.sim_mutation_rate * args.sim_locus_length
            recom_0 = 4.0 * args.sim_popsize_Hum * args.sim_recombination_rate * args.sim_locus_length
            simulation_command = '%d 3 -t %.8f -r %.8f %d -I 3 1 1 1 0 ' %(args.sim_replicates, theta_0, recom_0, args.sim_locus_length)
            for i in range(4):
                simulation_command += '-n %d %.8f ' %(i + 1, [args.sim_popsize_Hum])"""
        elif args.sim_model == 'kuhlwilm2016':
            #Using a mutation rate of 0.5e-9 per year, and the aggregated results (Table S14)
            pass
            """
            mutation_rate = 0.5e-9 #mutations per base per year
            recombination_rate = 1.25e-8 #recombination per base per year... check if right
            generation_time = 29.0
            N_0 = 1000.0 #arbitrary, as in Kuhlwilm
            Ne_H, Ne_N_70_140, Ne_N_140_450, Ne_D_50_450, Ne_DN, Ne_DNH  = [24000.0, 700.0, 3400.0, 2500.0, 7100.0, 18500.0]
            T_DN, T_DNH = [411700.0 / generation_time, 555200.0 / generation_time]
            S_D, S_N = [50000.0 / generation_time, 70000.0 / generation_time]
            theta_0 = 4.0 * N_0 * mutation_rate * args.sim_locus_length
            recom_0 = 4.0 * N_0 * recombination_rate * args.sim_locus_length
            introgression = 2.0
            simulation_command = '%d 3 -T -seed %f -t %f -r %f %d -I 3 1 1 1 0 -n 1 %f -n 2 %f -n 3 %f' %(time.time(),
                                                                                                          args.sim_replicates,
                                                                                                          theta_0,
                                                                                                          recom_0,
                                                                                                          args.sim_locus_length,
                                                                                                          Ne_H / float(N_0),
                                                                                                          Ne_N_70_140 / float(N_0),
                                                                                                          Ne_D_50_450 / float(N_0))
            simulation_command += ' -en %f %f -em %f 1 2 %f -em %f 1 2 0 -ej %f 3 2 -en %f 2 %f -ej %f 2 1 -en %f 1 %f' %((140000.0 / generation_time) / (4.0 * N_0),
                                                                                                                          Ne_N_140_450 / float(N_0),
                                                                                                                          migration time,
                                                                                                                          migration amount,
                                                                                                                          migration time - 1 generation,
                                                                                                                          T_DN / (4.0 * N_0),
                                                                                                                          T_DN / (4.0 * N_0) + 0.000001,
                                                                                                                          Ne_DN / float(N_0),
                                                                                                                          T_DNH / (4.0 * N_0),
                                                                                                                          T_DNH / (4.0 * N_0) + 0.000001,
                                                                                                                          Ne_DNH . float(N_0))
            #Need to include migration - potentially of both archaics populations. And specify that I sample before the present.

                                                                            
        """
        elif args.sim_model == 'malaspinas2016':
            generation_time = 29.0
            mutation_rate = 1.25e-8  #mutations per base per 29 year generation
            recombination_rate = 1.12e-8 #recombination per base per 29 year generation
            
            N_0 = 1000.0
            simulation_command = "scrm 3 1 -seed %f -t %f -r %f %d -I 5 1 1 1 0 0 0.0 -n 1 %f -n 2 %f -n 3 %f -n 4 %f -n 5 %f -eps %f 1 4 %f -eps %f 1 5 %f -ej %f 4 2 -en %f 2 %f -ej %f 5 3 -en %f 2 %f -ej %f 3 2 -en %f 2 %f -ej %f 2 1 -en %f 2 %f -T" %(time.time(), mutation_rate * args.sim_locus_length * 4 * N_0, recombination_rate * args.sim_locus_length * 4 * N_0, args.sim_locus_length, 23275.0 / N_0, 2846.0 / N_0, 463.0 / N_0, 7419.0/N_0, 7419.0 / N_0, 1517.0 / (4 * N_0), 1.0 - 0.04, 2069.0 / (4 * N_0), 1.0 - 0.031, 13580.0 / (4 * N_0), 13580.5 / (4 * N_0), 18296.0/N_0, 3820.0 / (4 * N_0), 3820.5 / (4 * N_0), 18296.0/ N_0, 17080.0 / (4 * N_0), 17080.5 / (4 * N_0), 18296.0/ N_0, 22652.0 / (4 * N_0), 22652.5 / (4 * N_0), 18296.0/ N_0)
            print simulation_command
        #run the simulation using scrm
        sim_run_name = "power_simulation.%s.%s.%d.txt" %(args.sim_model,args.viterbi_method,int(time.time() * 1000))
        f = open("I:/Genetic_Data/Indonesia_Diversity175/multicall3_mask1/phase/VCFsWithArchaic/Allele_sharing/simtmp/%s" %(sim_run_name), 'w')
        sim_run = subprocess.Popen(args = simulation_command,
                                   executable = 'I:/PhylogenPrograms/scrm/scrm.exe',
                                   stdout = f,
                                   stderr = f)
        sim_run.communicate()
        f.close()
        #convert the Newick trees in each simulation to Deni, Nean or non-introgressed. In each case, also write the properties.
        
        raise RuntimeError()
        
    #For each target, for each chromosome, do archaic annotation, summarise the output, and then do HMM Viterbi method
    #Determine the full file list for calculating archaich chunks
    input_file_list = []
    if len(args.infile) > 1 or (len(args.infile) == 1 and len(args.chromosomes) == 1): #Multiple infile specified, or both a single infile and chromosome number specified
        input_file_list = args.infile
        chromosome_list = [int(chrom) for chrom in args.chromosomes]
        input_recom_list = args.infile_recom
    else: #Multiple chromosomes and one infile
        if len(args.chromosomes) == 3:
            if args.chromosomes[1] == '...':
                chromosome_list = range(int(args.chromosomes[0]), int(args.chromosomes[2]) + 1)
            else:
                chromosome_list = [int(chrom) for chrom in args.chromosomes]
        else:
            chromosome_list = [int(chrom) for chrom in args.chromosomes]
        input_file_list = [args.infile[0] %(chrom) for chrom in chromosome_list]
        input_recom_list = [] if len(args.infile_recom) == 0 else [args.infile_recom[0] %(chrom) for chrom in chromosome_list]
    #female_list = [] #X chromosome not implemented
    #if args.females_file != 'None':
    #    with open(args.females_file, 'r') as f_fem:
    #        for line in f_fem:
    #            for ind in line[0:-1].split('\t'):
    #                female_list.append(ind)
    
    ##Read in the out/ingroup population, if one is provided.
    if args.popfile is not '':
        pop_inds = readin_populationFileEOL(args.popfile)
    else:
        pop_inds = []
    
    ##Determine the classification scheme
    if args.viterbi_method == 'binary_archaic_derived_sharing':
        if (len(args.viterbi_init_parameters) != (3 * args.viterbi_init_alternatives) and args.viterbi_init_resamples == 1) or (args.viterbi_init_resamples > 1 and len(args.viterbi_init_parameters) != 6):
            raise RuntimeError("binary_archaic_derived_sharing expects three parameters in args.viterbi_init_parameters per initial condition, or 6 for sampled - the initial proportion of the genome expected to be archaic, the rate of archiac signal in introgressed regions, and the rate of archaic signal in non-introgressed regions")
        classification_scheme_A = [[[[4,0],[6,3]]],[[[4,1],[6,1]]]] #0 if different vs archaic (4,0) and target is ancestral (6,3); 1 if same as archaic (4,1) and target is derived (6,1)
        #X chromosome not implemented
        classification_scheme_B = [[[[5,0],[7,3]]],[[[5,1],[7,1]]]]
        exclusion_scheme = [[[3,2]],[[10,1]], [[10,-1]]] #The allele shouldn't be '2' in the ancestral (3,2) or shared between the archaics (10,1), or missing in one archaic (10,-1)
    elif args.viterbi_method == 'binary_archaic_derived_sharing_no_derived_divergence':
        if (len(args.viterbi_init_parameters) != (3 * args.viterbi_init_alternatives) and args.viterbi_init_resamples == 1) or (args.viterbi_init_resamples > 1 and len(args.viterbi_init_parameters) != 6):
            raise RuntimeError("binary_archaic_derived_sharing_no_derived_divergence expects three parameters in args.viterbi_init_parameters per initial condition, or 6 for sampled - the initial proportion of the genome expected to be archaic, the rate of archiac signal in introgressed regions, and the rate of archaic signal in non-introgressed regions")
        #The classification scheme is saying:
        #0 if not shared with archaic (4,0), target is ancestral (6,3) and the comparative archaic is different vs the target one (10,0) [0,1,0]
        #0 if not shared with archaic (4,0), target is derived (6,4) and the comparative archaic is the same as the target one (10,1) [1,0,0]
        #1 if shared with archaic (4,1), target is derived (6,1) and the comparative archaic is different vs the target one (10,0) [1,1,0]
        classification_scheme_A = [[[[4,0],[6,3],[10,0]],[[4,0],[6,4],[10,1]]],[[[4,1],[6,1],[10,0]]]] #0:Targ=0,Arc=1,Comp=0;Targ=1,Arc=0,Comp=0   1:Targ=1,Arc=1,Comp=0
        #X chromosome not implemented
        classification_scheme_B = [[[[5,0],[7,3],[10,0]],[[5,0],[7,4],[10,1]]],[[[5,1],[7,1],[10,0]]]]
        exclusion_scheme = [[[3,2]],[[10,-1]]] #The allele shouldn't be '2' in the ancestral, or missing in either. To get 0:Targ=1,Arc=0,Comp=0 we can't remove all sharing, which means sharing conditions above
    elif args.viterbi_method == 'binary_archaic_derived_sharing_no_derived_divergence_joint_divergence':
        if (len(args.viterbi_init_parameters) != (3 * args.viterbi_init_alternatives) and args.viterbi_init_resamples == 1) or (args.viterbi_init_resamples > 1 and len(args.viterbi_init_parameters) != 6):
            raise RuntimeError("binary_archaic_derived_sharing_no_derived_divergence_joint_divergence expects three parameters in args.viterbi_init_parameters per initial condition, or 6 for sampled - the initial proportion of the genome expected to be archaic, the rate of archiac signal in introgressed regions, and the rate of archaic signal in non-introgressed regions")
        #The classification scheme is saying:
        #0 if not shared with archaic (4,0) and the target is ancestral (6,3). The state in the comparative archaic doesn't matter. [0,1,0], [0,1,1]
        #0 if not shared with archaic (4,0), target is derived (6,4) and the comparative archaic is the same as the target one (10,1) [1,0,0]
        #1 if shared with archaic (4,1), target is derived (6,1) and the comparative archaic is different vs the target one (10,0) [1,1,0]
        classification_scheme_A = [[[[4,0],[6,3]],[[4,0],[6,4],[10,1]]],
                                   [[[4,1],[6,1],[10,0]]]] #0:Targ=0,Arc=1,Comp=0;Targ=1,Arc=0,Comp=0   1:Targ=1,Arc=1,Comp=0
        #X chromosome not implemented
        classification_scheme_B = [[[[5,0],[7,3]],[[5,0],[7,4],[10,1]]],
                                   [[[5,1],[7,1],[10,0]]]]
        exclusion_scheme = [[[3,2]],[[10,-1]]] #The allele shouldn't be '2' in the ancestral, or missing in either. To get 0:Targ=1,Arc=0,Comp=0 we can't remove all sharing, which means sharing conditions above
    elif args.viterbi_method == 'binary_archaic_derived_sharing_racimo2016':
        print "\nUsing the test method binary_archaic_derived_sharing_racimo2016 . This is modelled on the emissions of Racimo et al 2016 TBX15/WARS2 paper, but there are several differences overall. a) The learning; they optimised vs simulations, I learn on data by EM. b) Possibly the algorith; their Git is Fwd-Back but SI says Viterbi. c) The SNP subset; they considered all SNPs polymorphic in the target population, I condsider all that have at least 1 derived in archaic or target humna. Note that I could consider polymorphic over the whole set, but this leads to lots of [0,0,x,>0] as signals against introgression. Similarly a signal [0,0,x,x] when the SNP is polymorphic in the target's population as against introgression isn't obvious. d) The transition probabilities are r*m*(t-1) and r*(1-m)*(t-1) for them, while I keep them linear but don't constrain them.\n"
        if (len(args.viterbi_init_parameters) != (3 * args.viterbi_init_alternatives) and args.viterbi_init_resamples == 1) or (args.viterbi_init_resamples > 1 and len(args.viterbi_init_parameters) != 6):
            raise RuntimeError("binary_archaic_derived_sharing_racimo2016 expects three parameters in args.viterbi_init_parameters per initial condition, or 6 for sampled - the initial proportion of the genome expected to be archaic, the rate of archiac signal in introgressed regions, and the rate of archaic signal in non-introgressed regions")
        #Racimo are 1 if [1,1,x,<min] and 0 if [0,1,x,x], [1,0,x,x], [1,1,x,>min] OR [0,0,x,x] if the SNP is polymorphic in the 'target population' i.e. MPI if individual is MPI.
        #The final category doesn't really make sense - it is saying that variation in the *population as a whole* is evidence against introgression.
        #I'm more comfortable doing things this way:
        #0 if [0,1,x,x], [1,0,x,x], [1,1,x,>min]
        #1 if [1,1,x,<min]
        #i.e.
        #0 if different from archaic (4,0) or shared derived with archaic (6,1) and freq in Africa > eta (12, minf + 1e-9, 1.0)
        #1 if shared with archaic and archaic target is derived (6,1), and frequency in Africa is < eta (12,0.0, minf).
        minf = (1.0 / (len(pop_inds) * 2)) + 0.00001
        classification_scheme_A = [[[[4,0]], [[6,1],[12,minf+1e-9, 1.0]]],
                                   [[[6,1], [12,0.0,minf]]]]
        
        classification_scheme_B = [[[[5,0]], [[7,1],[12,minf+1e-9, 1.0]]],
                                   [[[7,1], [12,0.0,minf]]]]
        
        exclusion_scheme = [[[3,2]],[[3,'.']],[[3,'nan']]] #The allele shouldn't be '2' in the ancestral, or be '.' in ancestral.
    
    elif args.viterbi_method == 'six_category_simple_archaic_demography':
        #The classification scheme says, with [Hum, Archaic, OutArchaic]:
        #0 if [1,0,0] = not shared with archaic and target derived (6,4), and archiacs are the same (10,1)
        #1 if [0,1,0] = not shared with archaic and target ancestral (6,3), and archiacs are different (10,0)
        #2 if [0,0,1] = shared with archaic and target ancestral (6,0), and archiacs are different (10,0)
        #3 if [1,1,0] = shared with archaic and target derived (6,1), and archiacs are different (10,0)
        #4 if [1,0,1] = not shared with archaic and target derived (6,4), and archiacs are different (10,0)
        #5 if [0,1,1] = not shared with archaic and target ancestral (6,3), and archiacs are the same (10,1)
        classification_scheme_A = [[[[6,4], [10,1]]], [[[6,3],[10,0]]], [[[6,0],[10,0]]], [[[6,1],[10,0]]], [[[6,4],[10,0]]], [[[6,3],[10,1]]]]
        classification_scheme_B = [[[[7,4], [10,1]]], [[[7,3],[10,0]]], [[[7,0],[10,0]]], [[[7,1],[10,0]]], [[[7,4],[10,0]]], [[[7,3],[10,1]]]]
        exclusion_scheme = [[[3,2]],[[10,-1]]]#,[[2,0],[3,1]],[[2,1],[3,0]]] #The allele shouldn't be '2' in the ancestral, or missing in either. Technically may be able to include '2' in one ancestral somehow, but leave for now.
    elif args.viterbi_method == 'fourteen_category_explore_archaic_pop_sharing':
        #The classification scheme is, with [Hum, Archaic, OutArchaic, OutPop]:
        #0 if [1,0,0,0], 1 if [0,1,0,0], 2 if [0,0,1,0], 3 if [0,0,0,1]
        #4 if [1,1,0,0], 5 if [1,0,1,0], 6 if [1,0,0,1], 7 if [0,1,1,0], 8 if [0,1,0,1], 9 if [0,0,1,1]
        #10 if [1,1,1,0], 11 if [1,1,0,1], 12 if [1,0,1,1], 13 if [0,1,1,1]
        minf = 1e-6
        classification_scheme_A = [[[[6,4],[10,1],[12,-0.1,0.0]]], [[[6,3],[10,0],[12,-0.1,0.0]]], [[[6,0],[10,0],[12,-0.1,0.0]]], [[[6,0],[10,1],[12,minf,1.0]]],
                                   [[[6,1],[10,0],[12,-0.1,0.0]]], [[[6,4],[10,0],[12,-0.1,0.0]]], [[[6,4],[10,1],[12,minf,1.0]]], [[[6,3],[10,1],[12,-0.1,0.0]]], [[[6,3],[10,0],[12,minf,1.0]]], [[[6,0],[10,0],[12,minf,1.0]]],
                                   [[[6,1],[10,1],[12,-0.1,0.0]]], [[[6,1],[10,0],[12,minf,1.0]]], [[[6,4],[10,0],[12,minf,1.0]]], [[[6,3],[10,1],[12,minf,1.0]]]]
        classification_scheme_B = [[[[7,4],[10,1],[12,-0.1,0.0]]], [[[7,3],[10,0],[12,-0.1,0.0]]], [[[7,0],[10,0],[12,-0.1,0.0]]], [[[7,0],[10,1],[12,minf,1.0]]],
                                   [[[7,1],[10,0],[12,-0.1,0.0]]], [[[7,4],[10,0],[12,-0.1,0.0]]], [[[7,4],[10,1],[12,minf,1.0]]], [[[7,3],[10,1],[12,-0.1,0.0]]], [[[7,3],[10,0],[12,minf,1.0]]], [[[7,0],[10,0],[12,minf,1.0]]],
                                   [[[7,1],[10,1],[12,-0.1,0.0]]], [[[7,1],[10,0],[12,minf,1.0]]], [[[7,4],[10,0],[12,minf,1.0]]], [[[7,3],[10,1],[12,minf,1.0]]]]
        exclusion_scheme = [[[3,2]],[[10,-1]]]
    else:
        raise RuntimeError("%s is not an implemented method" %(args.viterbi_method))
    
    #I now permit running the HMM multiple times using different initial conditions.
    #The initial conditions are generated early, and applied identically to each individual.
    if args.viterbi_init_resamples > 1 and args.viterbi_init_alternatives == 1:
        #In this case, args.viterbi_init_parameters is double the length of the number of intitial conditions needed, with each pair indicating the lower and upper bounds of a uniform distribution
        initial_conditions = []
        if len(args.viterbi_init_parameters) % 2 != 0:
            raise RuntimeError("Length of viterbi_init_parameters must be double the number of initial conditions if viterbi_init_resamples is used, current length is %d" %(len(args.viterbi_init_parameters)))
        for i in range(args.viterbi_init_resamples):
            initial_condition = []
            for j in range(int(len(args.viterbi_init_parameters) / 2)):
                lower_lim = args.viterbi_init_parameters[(j * 2)]
                upper_lim = args.viterbi_init_parameters[(j * 2) + 1]
                initial_condition.append((np.random.rand() * (upper_lim - lower_lim)) + lower_lim)
            initial_conditions.append(initial_condition)
    elif args.viterbi_init_alternatives > 1 and args.viterbi_init_resamples == 1:
        initial_conditions = []
        if len(args.viterbi_init_parameters) % args.viterbi_init_alternatives != 0:
            raise RuntimeError("Length of viterbi_init_parameters should be a multiple of viterbi_init_alternatives, currently %d" %(len(args.viterbi_init_parameters)))
        num_parameters = len(args.viterbi_init_parameters) / args.viterbi_init_alternatives
        for i in range(args.viterbi_init_alternatives):
            initial_conditions.append([args.viterbi_init_parameters[i + (j * args.viterbi_init_alternatives)] for j in range(num_parameters)])
    else:
        initial_conditions = [args.viterbi_init_parameters]
    print "Intial conditions including re-sample: ", initial_conditions
    
    for individual in range(len(args.target_individuals)):
        print "Estimating archaic haplotypes for individual %s" %(args.target_individuals[individual])
        print "Exclusion scheme is: ", exclusion_scheme
        print "Classification scheme (Chrom A) is: ", classification_scheme_A
        print "Classification scheme (Chrom B) is: ", classification_scheme_B
        pop_outfile = args.outfile if len(args.target_individuals) == 1 else args.outfile %(args.target_individuals[individual])
        log_outfile = args.outfile + '_log' if len(args.target_individuals) == 1 else args.outfile %(args.target_individuals[individual]) + '_log'
        profile_outfile = pop_outfile + '_profile'
        category_profile = [[],[]] #only used if profilemode is True
        gmap_file_list = []
        if args.log == True:
            fileoperation_log_append_from_list(log_outfile, line_list = [[],[],['\t'.join(["Individual:", str(args.target_individuals[individual])]),
                                                                                '\t'.join(["Comparison:", str(args.target_archaic)]),
                                                                                '\t'.join(["Viterbi method:", str(args.target_archaic)]),
                                                                                '\t'.join(["Viterbi transition parameter(s):"] + ['%.12f' %(i) for i in np.array([args.viterbi_transition]).flat]),
                                                                                '\t'.join(["Viterbi initial parameter(s):"] +  ['%.12f' %(i) for i in np.array(initial_conditions).flat ]),
                                                                                '\t'.join(["Exclusion scheme is:", str(exclusion_scheme)]),
                                                                                '\t'.join(["Classification scheme (Chrom A) is:", str(classification_scheme_A)]),
                                                                                '\t'.join(["Classification scheme (Chrom B) is:", str(classification_scheme_B)])]])
        for chrom_calculation in range(len(chromosome_list)):
            #Archaic annotation for the individual for that chromosome
            pop_archaic_annotation_outfile_chr = pop_outfile + '_archaicAnnotation_%s' %(chromosome_list[chrom_calculation])
            pop_archaic_annotation_summary_outfile_chr = pop_outfile + '_archaicAnnotationSummary_%s' %(chromosome_list[chrom_calculation])
            fileoperation_retrieve_archaic_admixture_annotation_infoaa(infile = input_file_list[chrom_calculation],
                                                                outfile = pop_archaic_annotation_outfile_chr,
                                                                target_name = args.target_individuals[individual],
                                                                reference = args.target_archaic,
                                                                nean_name = args.neanderthal_name,
                                                                deni_name = args.denisovan_name,
                                                                #anc_name = args.ancestral_name,
                                                                pop_list = pop_inds)#['MTW007','MTW020','MTW030','MTW057','MTW058','MTW062','MTW071','MTW078'])
            
            
            for chrom in [0,1]:
                #Classify each of the chromosome copies
                print "Converting archaic annotation data to input class string for Viterbi estimation chromosome %s" %('%d' %(chromosome_list[chrom_calculation]) + '_%d' %(chrom + 1))
                snps_posstart_posend = fileoperation_classify_archaic_annotation_as_string_array(infile = pop_archaic_annotation_outfile_chr,
                                                                    outfile_base = pop_archaic_annotation_summary_outfile_chr + '_%d' %(chrom + 1),
                                                                    classification_scheme = [classification_scheme_A, classification_scheme_B][chrom],
                                                                    exclusion_scheme = exclusion_scheme,
                                                                    check_conflicts = True,
                                                                    freq_columns = [11,12])
                if args.profilemode is False:
                    #If a genetic map is given, read it in.
                    if len(input_recom_list) != 0:
                        gmap_out = pop_archaic_annotation_summary_outfile_chr + '_%d_gmap' %(chrom + 1)
                        gmap_file_list.append(gmap_out)
                        try:
                            #I tend to have the gmap in two formats, either pos ?? cM...
                            fileoperation_interpolate_gmap_1KG(mapfile = input_recom_list[chrom_calculation], posfile = pop_archaic_annotation_summary_outfile_chr + '_%d_pos' %(chrom + 1), outfile = gmap_out, g_map_pos_IDx = 0, g_map_map_IDx = 2, g_map_delimiter = ' ')
                        except IndexError:
                            #or chr pos ?? cM
                            fileoperation_interpolate_gmap_1KG(mapfile = input_recom_list[chrom_calculation], posfile = pop_archaic_annotation_summary_outfile_chr + '_%d_pos' %(chrom + 1), outfile = gmap_out, g_map_pos_IDx = 1, g_map_map_IDx = 3, g_map_delimiter = '\t')
                    else:
                        gmap_out = None
                    chrom_transition = np.reshape(args.viterbi_transition, (int(np.sqrt(len(args.viterbi_transition))), int(np.sqrt(len(args.viterbi_transition)))))
                    print "Transition probability = ", chrom_transition
                    #Now call the Viterbi method.
                    print "Performing Viterbi optimisation using method %s" %(args.viterbi_method)
                    
                    if args.viterbi_method == 'binary_archaic_derived_sharing' or args.viterbi_method == 'binary_archaic_derived_sharing_no_derived_divergence' or args.viterbi_method == 'binary_archaic_derived_sharing_no_derived_divergence_joint_divergence' or args.viterbi_method == 'binary_archaic_derived_sharing_racimo2016':
                        curr_viterbi = None
                        curr_prob = -np.infty
                        for conditions in range(len(initial_conditions)):
                            print initial_conditions[conditions], initial_conditions[conditions][0]
                            proposed_viterbi_out = calculate_viterbiHMM_binary_archaic_derived_sharing(infile_types = pop_archaic_annotation_summary_outfile_chr + '_%d_types' %(chrom + 1),
                                                                            init_prob_high = np.log(initial_conditions[conditions][0]),
                                                                            init_prob_low_archaic = np.log(initial_conditions[conditions][1]),
                                                                            init_prob_high_archaic = np.log(initial_conditions[conditions][2]),
                                                                            transition_probability = np.log(chrom_transition),
                                                                            iterations = 50,
                                                                            fixed_parameters = args.viterbi_fixed_parameters,
                                                                            fit_transitions = args.viterbi_fit_transition_rates,
                                                                            genetic_map = gmap_out) #gmap_out is None if no gmap is supplied
                            proposed_viterbi_prob = support_calculate_viterbi_probabililty(state_path = proposed_viterbi_out[0], sequence_of_observations = proposed_viterbi_out[1], prob_low_archaic = proposed_viterbi_out[2][0], prob_high_archaic = proposed_viterbi_out[2][1], transition_probability = np.reshape(proposed_viterbi_out[4], (int(np.sqrt(len(proposed_viterbi_out[4]))),int(np.sqrt(len(proposed_viterbi_out[4]))))))
                            print proposed_viterbi_prob
                            print proposed_viterbi_out[2],proposed_viterbi_out[2][0]
                            print "Probability for initial conditions %s = %.3f . P(1|low introgressive class) = %.3f P(1|high introgressive class) = %.3f" %(' '.join(['%.4f' %(i) for i in initial_conditions[conditions]]),
                                                                                                                                                              proposed_viterbi_prob,
                                                                                                                                                              np.exp(proposed_viterbi_out[2][0]),
                                                                                                                                                              np.exp(proposed_viterbi_out[2][1]))
                            if proposed_viterbi_prob > curr_prob:
                                curr_viterbi = copy.copy(proposed_viterbi_out)
                                curr_prob = copy.copy(proposed_viterbi_prob)
                        viterbi_out = copy.copy(curr_viterbi)
                    if args.viterbi_method == 'six_category_simple_archaic_demography':
                        raise RuntimeError("Implement!")
                    else:
                        pass
                    #Write out the output as a list of tracts of each type, and write a log of the estimated parameters.
                    f_out_vit = open(pop_archaic_annotation_summary_outfile_chr + '_%d_types_viterbi' %(chrom + 1), 'wb')
                    f_out_vit.write(''.join(['%d' %(i) for i in viterbi_out[0]]) + '\n')
                    f_out_vit.close()
                    hap_lengths = [[] for i in range(int(max(viterbi_out[0])) + 1)]
                    prev_val, curr_type, curr_len = -1, -1, -1
                    for val in viterbi_out[0]:
                        if prev_val == -1:
                            prev_val = val
                            curr_type = val
                            curr_len = 1
                        elif val == prev_val:
                            curr_len += 1
                        else:
                            hap_lengths[int(curr_type)].append(curr_len)
                            prev_val = val
                            curr_type = val
                            curr_len = 1
                    if args.log == True:
                        ind_chrom_name = args.target_individuals[individual] + '_' + '%d' %(chromosome_list[chrom_calculation]) + '_%d' %(chrom + 1)
                        fileoperation_log_append_from_list(log_outfile, line_list = [[],['\t'.join(["%s" %(ind_chrom_name), "Total SNPs:", str(len(viterbi_out[0]))]),
                                                                                                  '\t'.join(["%s" %(ind_chrom_name), "Class_0-N_SNPs"] + ["%d" %(i) for i in viterbi_out[3]]),
                                                                                                  '\t'.join(["%s" %(ind_chrom_name), "Inferred_archaic_SNPs", str(np.sum(viterbi_out[0]))]),
                                                                                                  '\t'.join(["%s" %(ind_chrom_name), "Total_switchpoints", str(np.sum(np.abs(viterbi_out[0] - np.roll(viterbi_out[0], 1))) - 1)]),
                                                                                                  '\t'.join(["%s" %(ind_chrom_name), "Transition_probability", '%s' %(' '.join(['%.3f' %(i) for i in chrom_transition.flatten()]))]),
                                                                                                  '\t'.join(["%s" %(ind_chrom_name), "Average_haplotype_length", '%s' %(' '.join(['%.1f' %(np.average(i)) for i in hap_lengths]))]),
                                                                                                  '\t'.join(["%s" %(ind_chrom_name), "Inferred_parameters", '%s' %(' '.join(["%.4f" %(i) for i in np.exp(viterbi_out[2])]))]),
                                                                                                  '\t'.join(["%s" %(ind_chrom_name), "Inferred_transitions", '%s' %(' '.join(["%.4f" %(i) for i in np.exp(viterbi_out[4])]))])]])
                    print "Inferred %d of %d SNPs in the high archaic class with %d switchpoints and average haplotype lengths [%s] median [%s] SNPs by type. Inferred parameters are [%s]. Inferred transitions are [%s]." %(np.sum(viterbi_out[0]),
                                                                                                                                                                                   len(viterbi_out[0]),
                                                                                                                                                                                   np.sum(np.abs(viterbi_out[0] - np.roll(viterbi_out[0], 1))) - 1,
                                                                                                                                                                                   ' '.join('%.1f' %(np.average(i)) for i in hap_lengths),
                                                                                                                                                                                   ' '.join('%.1f' %(np.median(i)) for i in hap_lengths),
                                                                                                                                                                                   ' '.join("%.6f" %(i) for i in np.exp(viterbi_out[2])),
                                                                                                                                                                                   ' '.join(["%.12f" %(i) for i in np.exp(viterbi_out[4])]))

                elif args.profilemode is True:
                    #Just count the types and save.
                    with open(pop_archaic_annotation_summary_outfile_chr + '_%d_types' %(chrom + 1), 'rb') as f:
                        line = f.readline()[0:-1]
                    categories = np.array([line.count('%s' %(chr(i))) for i in range(48,91)])
                    category_profile[chrom].append(categories)
                else:
                    pass
            if args.cleanup_level >= 1:
                os.remove(pop_archaic_annotation_outfile_chr)
        if args.profilemode is True:
            category_profile = np.array(category_profile)
            #print category_profile
            max_cat = max(np.where(np.sum(category_profile[0], 0) > 0)[0][-1], np.where(np.sum(category_profile[1], 0) > 0)[0][-1])
            for chrom in [0,1]:
                with open(profile_outfile + '_%d' %(chrom + 1), 'wb') as f:
                    for i in range(len(category_profile[chrom])):
                        line_to_write = ['%d' %(chromosome_list[i])] + ['%d' %(j) for j in category_profile[chrom][i][0:max_cat + 1]]
                        #print line_to_write
                        f.write('\t'.join(line_to_write) + '\n')
                    f.write('total\t' + '\t'.join(['%d' %(j) for j in np.sum(category_profile[chrom], 0)[0:max_cat + 1]]) + '\n')

        #End of chrom list for individual. Combine Viterbi into BED of 0-chunks and 1-chunks
        elif args.chrom_to_bed is True:
            for chrom in [0,1]:
                tmp = pop_outfile + '_archaicAnnotationSummary_%s' + '_%d' %(chrom + 1)
                fileoperation_viterbi_positions_to_BED(filebase = pop_outfile + '_archaicAnnotationSummary_%s' + '_%d' %(chrom + 1),
                                                       outfile = pop_outfile + '_archaicAnnotationSummary_combined_types_%d' %(chrom + 1),
                                                       type_suffix = '_types',
                                                       chroms = ['%d' %(i) for i in chromosome_list],
                                                       track_name = args.target_individuals[individual] + '_%d' %(chrom + 1),
                                                       track_description = 'Exclusion:%s Classification:%s' %(exclusion_scheme, classification_scheme_A))
                fileoperation_viterbi_positions_to_BED(filebase = pop_outfile + '_archaicAnnotationSummary_%s' + '_%d' %(chrom + 1),
                                                       outfile = pop_outfile + '_archaicAnnotationSummary_combined_typesViterbi_%d' %(chrom + 1),
                                                       type_suffix = '_types_viterbi',
                                                       chroms = ['%d' %(i) for i in chromosome_list],
                                                       track_name = args.target_individuals[individual] + '_%d' %(chrom + 1),
                                                       track_description = 'Exclusion:%s Classification:%s' %(exclusion_scheme, classification_scheme_B))
        if args.cleanup_level >= 2:
            for chrom in [0,1]:
                for chrom_calculation in range(len(chromosome_list)):
                    os.remove(pop_outfile + '_archaicAnnotationSummary_%d_%d_pos' %(chromosome_list[chrom_calculation], chrom + 1))
                    os.remove(pop_outfile + '_archaicAnnotationSummary_%d_%d_types' %(chromosome_list[chrom_calculation], chrom + 1))
                    if args.profilemode is False:
                        os.remove(pop_outfile + '_archaicAnnotationSummary_%d_%d_types_viterbi' %(chromosome_list[chrom_calculation], chrom + 1))
            for rm_gfile in gmap_file_list:
                os.remove(rm_gfile)
