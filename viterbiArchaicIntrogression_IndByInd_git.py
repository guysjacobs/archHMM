###This is a program to estimate Denisovan/Neanderthal ancestry haplotypes in a target
###individual.

###The implemetation was designed for flexibility, allowing the detection of tracts
###conditioned on multiple archaic hominins for example, or without conditioning on
###a human outgroup. However, for ease of use this version just includes the HMM
###described in Jacobs et al 2018.

###The sequence of operations is:
###a) annotate each SNP based on its allele in target individual, archaic, outgroup etc.
###b) classify each SNP based its annotation into a string of 1s and 0s
###c) call the Viterbi algorithm on this string, using EM to estimate transition and emission probs
###d) convert the output to a BED file, indicating the span of inferred haplotypes, where they start and end, but not the specific SNPs that contribute.

###It is quite slow (e.g. ~one hour per genome).

import sys
import numpy as np
import argparse
import os
import subprocess
import time
import gzip

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/functions/')
print sys.path
from fileoperation_classify_archaic_annotation_as_string_array import *
from fileoperation_retrieve_archaic_admixture_annotation_infoaa import *
from fileoperation_log_append_from_list import *
from fileoperation_interpolate_gmap_1KG import *
from calculate_viterbiHMM_binary_archaic_derived_sharing import calculate_viterbiHMM_binary_archaic_derived_sharing, support_calculate_viterbi_probabililty
from fileoperation_viterbi_positions_to_BED import *
from readin_populationFileEOL import *

VERSION = 0.21

parser = argparse.ArgumentParser(description='Calculate tracts of archaic introgression based on a Viterbi algorithm and allele sharing, based on a VCF file including both a target individual and archaics.')

parser.add_argument('--infile', metavar='infile', type=str, nargs='*', default = ['None'],
                    help='system location of the VCF input file; multiple e.g. chromosome files can be passed, but if so their respective chromosomes must still be specified using --chromosomes.')
parser.add_argument('--outfile', metavar='outfile', type=str, nargs='?', default = 'None',
                    help='system location of the output file root. This program generates two outputs: _tracts (indicating when different tracts begin and end) and _probs (indicating the fitted archaic frequencies)')

parser.add_argument('--infile_recom', metavar='infile_recom', type=str, nargs='*', default = [],
                    help='system location of a genetic map input file, if available. Multiple e.g. chromosome files can be passed, but if so their respective chromosomes must still be specified using --chromosomes.')

parser.add_argument('--target_individuals', '-targi', dest='target_individuals', type=str, action = 'store', nargs='*', default = [],
                    help='which individuals to calculate archaic admixture tracts for; this is the name in the VCF.')
parser.add_argument('--target_archaic', '-targa', dest='target_archaic', type=str, action = 'store', nargs='?', default = '', choices = ['Neanderthal', 'Denisovan'],
                    help='estimate introgressing Denisovan or Neanderthal haplotypes?')
parser.add_argument('--denisovan_name', '-deni', dest='denisovan_name', type=str, action = 'store', nargs='?', default = '',
                    help='name of the target Denisovan individual in the VCF.')
parser.add_argument('--neanderthal_name', '-nean', dest='neanderthal_name', type=str, action = 'store', nargs='?', default = '',
                    help='name of the target Neanderthal individual in the VCF.')
parser.add_argument('--chromosomes', '-c', dest='chromosomes', action = 'store', nargs = '*', default = [],
                    help="list (or range, '1 ... 22') of chromosomes to calculate archaic tracts on. These correspond in order to the files given as infile, or are substituted as integers into a single infile at the '%%d' position.")
parser.add_argument('--log', '-l', dest='log', action = 'store_true', default = False,
                    help="should log files be written?")
parser.add_argument('--popfile', '-popfile', dest='popfile', type=str, action = 'store', nargs = '?', default = '',
                    help="population file (one individual name per line) to be used when allele classes are conditioned on population frequency")

parser.add_argument('--viterbi_method', '-vm', dest='viterbi_method', action = 'store', nargs='?', default = 'binary_archaic_derived_sharing_racimo2016', choices = ['binary_archaic_derived_sharing_racimo2016', 'manual_specification'],
                    help='the annotation/HMM method to be applied. The default operational option is currently binary_archaic_derived_sharing_racimo2016. Use manual_specification to read in a custom method.')
parser.add_argument('--viterbi_specification_file', '-vspecf', dest='viterbi_specification_file', action = 'store', nargs='?', default = '', type = str,
                    help='the manual HMM method specification file. This is only required if --viterbi_method is manual_specificaion. This file can be used to specify how the HMM runs - e.g. the minimum frequency, what SNP motifs count as what emssions 0 and 1. This method is primarily for testing different HMM approaches, adding some flexibilitiy...')


parser.add_argument('--viterbi_init_parameters', '-vparam', dest='viterbi_init_parameters', type=float, action = 'store', nargs='*', default = [],
                    help='the initial parameters for the viterbi algorithm. Here, we expect three parameters - the initial probability of a SNP being introgressed, the emission probability of archaic signal (1s) in low-archaic regions and the rate of archaic signal in high-archaic regions.')
parser.add_argument('--viterbi_fixed_parameters', '-vparamfix', dest='viterbi_fixed_parameters', type=int, action = 'store', nargs='*', default = [],
                    help='the IDxs of the parameters to keep fixed when optimising the chunks. This can be relevant for between-population comparisons or if a single introgression event is assumed. Be aware that this radically changes what the program is telling us. In the case of binary_archaic_derived_sharing and similar, the parameter IDxs are 0 for the initial amount of introgressive SNPs in the low class and 1 for the initial amount of introgressive SNPs in the high class.')

parser.add_argument('--viterbi_transition', '-vt', dest='viterbi_transition', type=float, action = 'store', nargs='*', default = [0.9,0.1,0.1,0.9],
                    help='the transition probability matrix for switching states. Low tranistion probabilities lead to longer haplotypes. The full matrix should be specified i.e. for two hidden states: a, b, c, d is [[a,b],[c,d]] is [[fromAtoA,fromAtoB], [fromBtoA,fromBtoB]]. a and d are the probabilities of no transition, and a should be 1.0 - b and d is 1.0 - c. I now offer an option to fit this parameter using viterbi_fit_transition_rates .')
parser.add_argument('--viterbi_fit_transition_rates', '-vfittrans', dest='viterbi_fit_transition_rates', action = 'store_true', default = False,
                    help="fit the transition rates as well as the probabilities. This involves counting transitions along the chromosome after each Viterbi iteration, and updating. Best used when simply searching for the best model for a specific chromosome. May have high sensitivity to initial conditions.")
parser.add_argument('--viterbi_EM_mode', '-vEMmode', dest='viterbi_EM_mode', action = 'store', type = str, nargs = '?', default = 'independent_chroms', choices = ['independent_chroms', 'save_fitted_parameters', 'readin_fitted_parameters'],
                    help='How to do the EM fitting? There are three options - each chromosome may be fitted independently (independent_chroms); you can do a run in which I specifically save the fitted parameters for all chromosomes included in the analysis (save_fitted_parameters); or you can do a run in which previously saved fitted parameters are read in and averaged, with these parameters then taken as fixed (readin_fitted_parameters). The latter two are used together by running save_fitted_parameters and then readin_fitted_parameters.')

parser.add_argument('--cleanup_level', '-cl', dest = 'cleanup_level', type = int, action = 'store', nargs = '?', default = '0',
                    help="which files to delete during the calculation. 0 means no files are deleted; 1 deletes the archaicAnnotation file; 2 additionally deletes the chromosome-by-chromosome output files; 3 additionally deletes any fitted emissions or transmissions files; 4 additionally deletes all results files (just left with logs).")
parser.add_argument('--chrom_to_bed', '-chrombed', dest='chrom_to_bed', action = 'store_true', default = True,
                    help="convert the Viterbi output files to a BED file for comparing individuals/ease of manipulation.")

args = parser.parse_args()

print 'Input file(s): ',  ' '.join(['%s' %(infile) for infile in args.infile])
print 'Output file: ', args.outfile
print 'Target individual list: ', ' '.join(['%s' %(ind) for ind in args.target_individuals])
print 'Chromosomes used: ', ' '.join(['%s' %(chrom) for chrom in args.chromosomes])
print "Initial conditions: ", args.viterbi_init_parameters
print "Transition inital conditions: ", args.viterbi_transition

if args.infile[0] == "None" or args.outfile == "None":
    print "Please specify input and output files."

elif len(args.target_individuals) == 0 or (args.target_archaic != 'Neanderthal' and args.target_archaic != 'Denisovan'):
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

elif len(args.infile_recom) == 0 and np.sum(np.sum(np.reshape(args.viterbi_transition, (int(np.sqrt(len(args.viterbi_transition))), int(np.sqrt(len(args.viterbi_transition))))), 1) == np.ones(int(np.sqrt(len(args.viterbi_transition))))) != int(np.sqrt(len(args.viterbi_transition))):
    print "viterbi_transition is the transition matrix, which are probabilities if no recombination map is supplied and so should sum to K (number of hidden states)."

elif len(args.viterbi_transition) != 4:
    print "viterbi_transition is the transition matrix. As there are two hidden states it must have 4 entries."

elif args.viterbi_method == 'manual_specification' and args.viterbi_specification_file == '':
    print "viterbi_method is manual_specification. This requires an emission specification file which is provided using --viterbi_specification_file PATH_TO_SPECIFICATION_FILE"

#TESTING I:\Dropbox\PythonPrograms\github\archHMM>python ./viterbiArchaicIntrogression_IndByInd_git.py --infile I:\Dropbox\Transfer\Indonesia_Diversity175_data\VCFsNew\chr%d.AD.AN_highQ_only_wo_N_biSNP.ANC.479.haps.vcf.gz --outfile I:\Dropbox\Transfer\Indonesia_Diversity175_data\hmm\herrInDeni\herrInDeni_CaseA_%s --infile_recom I:\Dropbox\Transfer\Indonesia_Diversity175_data\Genetic_maps\genetic_map_HapMapII_GRCh37\genetic_map_GRCh37_chr%d.txt.gz --chromosomes 21 22 --neanderthal_name AltaiNea --denisovan_name DenisovaPinky --target_individuals DenisovaPinky --target_archaic Neanderthal --log --popfile I:\Dropbox\Transfer\Indonesia_Diversity175_data\population_data\sample_lists\list_subset_africa_SSonly.txt --viterbi_method manual_specification --viterbi_specification_file I:\Dropbox\PythonPrograms\github\manual_specification_HMMHErDeni_CaseA.txt --viterbi_init_parameters 0.1 0.05 0.2 --viterbi_transition 0.9 0.1 0.1 0.9 --viterbi_fit_transition_rates --viterbi_EM_mode independent_chroms --cleanup_level 0

else:
    #For each target, for each chromosome, do archaic annotation, summarise the output, and then do HMM Viterbi method
    #Determine the full file list for calculating archaich chunks
    if args.cleanup_level >= 4:
        print "WARNING: This mode generates results, but then deletes them. This is intended for use when args.viterbi_EM_mode == save_fitted_parameters, as the results of these runs may not be used"
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
    
    ##Read in the out/ingroup population, if one is provided.
    if args.popfile is not '':
        pop_inds = readin_populationFileEOL(args.popfile)
    else:
        pop_inds = []
    
    ##Determine the classification scheme.
    ##For simplicity, the only classification scheme option is currently binary_archaic_derived_sharing_racimo2016
    if args.viterbi_method == 'binary_archaic_derived_sharing_racimo2016':
        print "\nUsing the method binary_archaic_derived_sharing_racimo2016 . This is similar to the Racimo et al 2016 TBX15/WARS2 paper, but there are several differences. a) The learning; this program tries to learn on data by EM. b) The algorithm is Viterbi, the Racimo 2016 Git has Fwd-Back. c) The SNP subset; I think that they considered all SNPs that are polymorphic in the target population, I only consider SNPs all that have at least 1 derived in archaic or target human. d) The transition probabilities are r*m*(t-1) and r*(1-m)*(t-1) for them, while I keep them linear but don't constrain them.\n"
        if len(args.viterbi_init_parameters) != 3:
            raise RuntimeError("binary_archaic_derived_sharing_racimo2016 expects three parameters in args.viterbi_init_parameters. The three paramters are the initial proportion of the genome expected to be archaic, the emission rate of archaic signal (1s) in less archaic regions, and the emission rate of archaic signal (1s) in more archaic regions")
        #My understanding of Racimo et al 2016 has a 1 for motif 1,1,x,<min] and 0 for motifs [0,1,x,x], [1,0,x,x], [1,1,x,>min] OR [0,0,x,x], for all SNP that are polymorphic in the 'target population'.
        #If that is correct then the final category implies that high variation in the *population as a whole* can be evidence against introgression.
        #Either way, I'm more comfortable doing things this way:
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
    elif args.viterbi_method == 'manual_specification':
        print "\nUsing a manually specified emission and exclusion scheme. This describes which mutation motifs correspond to which emissions and which SNPs should be excluded."
        classification_scheme_A = []
        classification_scheme_B = []
        exclusion_scheme = []
        insert_values = {}
        lambda_haps = [len(pop_inds) * 2.0]
        in_function = False
        f_exclusion = open(args.viterbi_specification_file, 'rb')
        try:
            for line in f_exclusion:
                if line[0:2] == '##':
                    pass
                elif line[0] == '#' or line[0] == '\n':
                    #new function
                    if in_function == True:
                        function_string = ''.join(function)
                        for key in insert_values.keys():
                            function_string = function_string.replace(key, str(insert_values[key]))
                        if function_name == 'classification_scheme_A':
                            classification_scheme_A = eval(function_string)
                        elif function_name == 'classification_scheme_B':
                            classification_scheme_B = eval(function_string)
                        elif function_name == 'exclusion_scheme':
                            exclusion_scheme = eval(function_string)
                        else:
                            lambda_function = eval('lambda x1 : ' + ''.join(function))
                            insert_values[function_name] = lambda_function(*lambda_haps)
                            
                    function_name = line[1:-1]
                    function = []
                    in_function = True
                else:
                    if in_function == True:
                        line_text = line
                        line_text = line_text.replace(' ', '')
                        line_text = line_text.replace('\n', '')
                        function.append(line_text)
        except:
            raise RuntimeError("Failed in manual setting of HMM emission and exlusion criteria. Check file formatting?")
    else:
        raise RuntimeError("%s is not an implemented method" %(args.viterbi_method))
    
    #For simplicity, this version only allows one set of initial conditions here at a time (rather than profiling multiple initial conditions).
    
    for individual in range(len(args.target_individuals)):
        print "Estimating archaic haplotypes for individual %s" %(args.target_individuals[individual])
        print "Exclusion scheme is: ", exclusion_scheme
        print "Classification scheme (Chrom A) is: ", classification_scheme_A
        print "Classification scheme (Chrom B) is: ", classification_scheme_B
        pop_outfile = args.outfile if len(args.target_individuals) == 1 else args.outfile %(args.target_individuals[individual])
        log_outfile = args.outfile + '_log' if len(args.target_individuals) == 1 else args.outfile %(args.target_individuals[individual]) + '_log'
        
        if args.viterbi_EM_mode == 'readin_fitted_parameters':
            print "Reading in fitted parameters. Note that this disables commands that would lead to repeating the EM fitting, i.e. the (previously fitted) parameters are taken as fixed"
            print "Setting args.viterbi_fit_transition_rates from %s to False" %(args.viterbi_fit_transition_rates)
            args.viterbi_fit_transition_rates = False
            print "Setting args.viterbi_fixed_parameters from %s to 0,1" %(','.join(str(i) for i in args.viterbi_fixed_parameters))
            args.viterbi_fixed_parameters = [0,1]
            with open(pop_outfile + '_fitted_emissions', 'rb') as f_emissions:
                fitted_emissions = []
                for line in f_emissions:
                    fitted_emissions.append(np.array(line[0:-1].split('\t')[1:], dtype = float))
                fitted_emissions = np.array(fitted_emissions)
                initial_conditions = [[args.viterbi_init_parameters[0], np.average(fitted_emissions[::,0]), np.average(fitted_emissions[::,1])]]
            with open(pop_outfile + '_fitted_transitions', 'rb') as f_transitions:
                fitted_transitions = []
                for line in f_transitions:
                    fitted_transitions.append(np.array(line[0:-1].split('\t')[1:], dtype = float))
                fitted_transitions = np.array(fitted_transitions)
                initial_transitions = [np.average(fitted_transitions[::,i]) for i in range(len(fitted_transitions[0]))]
        else:
            initial_conditions = [args.viterbi_init_parameters]
            initial_transitions = args.viterbi_transition
        print "Initial conditions: ", initial_conditions
        print "Initial transitions: ", initial_transitions
        estimated_chrom_emissions = []
        estimated_chrom_transitions = []
        
        gmap_file_list = []
        if args.log == True:
            fileoperation_log_append_from_list(log_outfile, line_list = [[],[],['\t'.join(["ScriptVersion:", str(os.path.basename(__file__))]),
                                                                                '\t'.join(["Input file(s):"] + args.infile),
                                                                                '\t'.join(["Individual:", str(args.target_individuals[individual])]),
                                                                                '\t'.join(["Comparison:", str(args.target_archaic), str(args.neanderthal_name) if args.target_archaic == 'Neanderthal' else str(args.denisovan_name) if args.target_archaic == 'Denisovan' else 'Misspecified']),
                                                                                '\t'.join(["Viterbi method:", str(args.viterbi_method)]),
                                                                                '\t'.join(["Viterbi transition parameter(s):"] + ['%.12f' %(i) for i in np.array([initial_transitions]).flat]),
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
                                                                pop_list = pop_inds)
            
            
            for chrom in [0,1]:
                #Classify each of the chromosome copies
                print "Converting archaic annotation data to input class string for Viterbi estimation chromosome %s" %('%d' %(chromosome_list[chrom_calculation]) + '_%d' %(chrom + 1))
                snps_posstart_posend = fileoperation_classify_archaic_annotation_as_string_array(infile = pop_archaic_annotation_outfile_chr,
                                                                    outfile_base = pop_archaic_annotation_summary_outfile_chr + '_%d' %(chrom + 1),
                                                                    classification_scheme = [classification_scheme_A, classification_scheme_B][chrom],
                                                                    exclusion_scheme = exclusion_scheme,
                                                                    check_conflicts = True,
                                                                    freq_columns = [11,12])
                #If a genetic map is given, read it in.
                if len(input_recom_list) != 0:
                    gmap_out = pop_archaic_annotation_summary_outfile_chr + '_%d_gmap' %(chrom + 1)
                    gmap_file_list.append(gmap_out)
                    try:
                        #I tend to have the gmap in two formats, either [pos,-,cM] with ' ' delimiter...
                        fileoperation_interpolate_gmap_1KG(mapfile = input_recom_list[chrom_calculation], posfile = pop_archaic_annotation_summary_outfile_chr + '_%d_pos' %(chrom + 1), outfile = gmap_out, g_map_pos_IDx = 0, g_map_map_IDx = 2, g_map_delimiter = ' ')
                    except IndexError:
                        #or [chr,pos,-,cM] with '\t' delimiter
                        try:
                            fileoperation_interpolate_gmap_1KG(mapfile = input_recom_list[chrom_calculation], posfile = pop_archaic_annotation_summary_outfile_chr + '_%d_pos' %(chrom + 1), outfile = gmap_out, g_map_pos_IDx = 1, g_map_map_IDx = 3, g_map_delimiter = '\t')
                        except:
                            raise RuntimeError("Check format of the genetic map. I expect either columns [pos,-,cM] with ' ' delimiter or [chr,pos,-,cM] with '\t' delimiter")
                else:
                    gmap_out = None
                chrom_transition = np.reshape(initial_transitions, (int(np.sqrt(len(initial_transitions))), int(np.sqrt(len(initial_transitions)))))
                print "Transition probability = ", chrom_transition
                
                #Now call the Viterbi method.
                print "Performing Viterbi optimisation using method %s" %(args.viterbi_method)
                
                if args.viterbi_method in ['binary_archaic_derived_sharing_racimo2016', 'manual_specification']:
                    curr_viterbi = None
                    curr_prob = -np.infty
                    for conditions in range(len(initial_conditions)):
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
                        print "Probability for initial conditions %s = %.3f . P(1|low introgression hidden state) = %.3f P(1|high introgression hidden state) = %.3f" %(' '.join(['%.4f' %(i) for i in initial_conditions[conditions]]),
                                                                                                                                                          proposed_viterbi_prob,
                                                                                                                                                          np.exp(proposed_viterbi_out[2][0]),
                                                                                                                                                          np.exp(proposed_viterbi_out[2][1]))
                        if proposed_viterbi_prob > curr_prob:
                            curr_viterbi = copy.copy(proposed_viterbi_out)
                            curr_prob = copy.copy(proposed_viterbi_prob)
                    viterbi_out = copy.copy(curr_viterbi)
                    estimated_chrom_emissions.append(['%d_%d' %(chromosome_list[chrom_calculation], chrom + 1), np.exp(proposed_viterbi_out[2][0]), np.exp(proposed_viterbi_out[2][1])])
                    estimated_chrom_transitions.append(['%d_%d' %(chromosome_list[chrom_calculation], chrom + 1), np.exp(proposed_viterbi_out[4][0]), np.exp(proposed_viterbi_out[4][1]),np.exp(proposed_viterbi_out[4][2]), np.exp(proposed_viterbi_out[4][3])])
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

            if args.cleanup_level >= 1:
                os.remove(pop_archaic_annotation_outfile_chr)
        
        #End of chrom list for individual. Combine Viterbi into BED of 0-chunks and 1-chunks
        if args.chrom_to_bed is True:
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
            #Remove chromosome-specific files
            for chrom in [0,1]:
                for chrom_calculation in range(len(chromosome_list)):
                    os.remove(pop_outfile + '_archaicAnnotationSummary_%d_%d_pos' %(chromosome_list[chrom_calculation], chrom + 1))
                    os.remove(pop_outfile + '_archaicAnnotationSummary_%d_%d_types' %(chromosome_list[chrom_calculation], chrom + 1))
                    os.remove(pop_outfile + '_archaicAnnotationSummary_%d_%d_types_viterbi' %(chromosome_list[chrom_calculation], chrom + 1))
            for rm_gfile in gmap_file_list:
                os.remove(rm_gfile)
        if args.cleanup_level >= 3:
            #Remove individual-specific fitted parameters files
            if os.path.isfile(pop_outfile + '_fitted_emissions'):
                os.remove(pop_outfile + '_fitted_emissions')
            else:
                pass
            if os.path.isfile(pop_outfile + '_fitted_transitions'):
                os.remove(pop_outfile + '_fitted_transitions')
            else:
                pass
        if args.cleanup_level >= 4:
            #Remove individual-specific combined results files!
            for chrom in [0,1]:
                os.remove(pop_outfile + '_archaicAnnotationSummary_combined_types_%d' %(chrom + 1))
                if os.path.isfile(pop_outfile + '_archaicAnnotationSummary_combined_typesViterbi_%d' %(chrom + 1)):
                    os.remove(pop_outfile + '_archaicAnnotationSummary_combined_typesViterbi_%d' %(chrom + 1))
                else:
                    pass

        #Write out the EM fitted parameters for all chromosomes included if args.viterbi_EM_mode == 'save_fitted_parameters'
        if args.viterbi_EM_mode == 'save_fitted_parameters':
            with open(pop_outfile + '_fitted_emissions', 'wb') as f_emissions:
                for chrom_hap in estimated_chrom_emissions:
                    f_emissions.write('%s\t%.12f\t%.12f\n' %(chrom_hap[0], chrom_hap[1], chrom_hap[2]))
            with open(pop_outfile + '_fitted_transitions', 'wb') as f_transitions:
                for chrom_hap in estimated_chrom_transitions:
                    f_transitions.write('%s\t%.12f\t%.12f\t%.12f\t%.12f\n' %(chrom_hap[0], chrom_hap[1], chrom_hap[2], chrom_hap[3], chrom_hap[4]))
        
