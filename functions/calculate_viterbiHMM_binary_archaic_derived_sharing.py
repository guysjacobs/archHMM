###This is an implementation of a HMM for detecting archaic admixture
###The method tries to categorise the genome into blocks of two hidden types, introgressed or not introgressed, given:
###A binary state vector of observed states
###Fixed or EM fitted emission probabilities
###Fixed or EM fitted transition probabilities

###The method can use a genetic map.

import numpy as np
import copy

def calculate_viterbiHMM_binary_archaic_derived_sharing(infile_types, init_prob_high, init_prob_low_archaic, init_prob_high_archaic, transition_probability, iterations = 100, fixed_parameters = [], fit_transitions = False, genetic_map = None):
    """
    This function calculates the viterbiHMM until it converges.
    
    The initial parameters are:
    init_prob_high, the initial probability of observing the introgressed hidden state
    init_prob_low_archaic, the initial probability of emmitting 1 given non-introgredded hidden state
    init_prob_high_archaic, the initial probability of emmitting 1 given introgredded hidden state

    The operations are:
    a) Read in the input string file, which has annotations for each (included) SNP as 0 (less likely given introgression) or 1 (more likely given introgression).
    b) Calculate the Viterbi algorithm using the initial emission and transmission probabilities.
    If requested (i.e. not all parameters are fixed, iterations > 1, and/or fit_transitions is True), do EM until convergence or we hit the iteration limit (iterations):
    c) Use the output to re-calculate the probability of each observations given each hidden state, and transitions.
    d) If the paramters have converged, stop; else feed the newly estimated parameters back in, re-calculate the Viterbi algorithm, and go to c)

    The process can be repeated multiple times with different parameters by inputting len>1 lists of values init_prob_high, init_prob_low_archaic, init_prob_high_archaic
    In this case, the iteration that has the highest log-probability is used.
    
    fixed_parameters specifies if any of the parameters are fixed at their initial value. In this case, 0 is for init_prob_low_archaic and 1 is for init_prob_high_archaic.
    
    There is an option to fit the transition probability using fit_transitions = True
    
    There is an option to incorporate a genetic map using genetic_map. The format should be one line per position, with two tab-delimited columns: position \t genetic map position
    """
    open_infile = gzip.open if infile_types[-3:] == '.gz' else open
    
    with open_infile(infile_types, 'rb') as f_in:
        #0-9 are characters and >9 are ascii to let me use up to 90-48 = 42 classes
        string_of_archaic_alleles = np.array(np.fromstring(f_in.readline()[0:-1], 'S1').view(np.ubyte) - 48, dtype = int)
    max_str = np.max([i if np.sum(string_of_archaic_alleles == i) > 0 else 0 for i in range(9)])
    print "Observed states are 0 to %d. Hidden states 0 and 1." %(max_str)
    
    if genetic_map is not None:
        gen_map_pos = []
        open_gmap = gzip.open if genetic_map[-3:] == '.gz' else open
        with open_gmap(genetic_map, 'rb') as f_gmap:
            for line in f_gmap:
                split_line = line[0:-1].split('\t')
                gen_map_pos.append(float(split_line[1])) #One line per SNP, genetic map position
        #Convert forward distance between SNPs
        gen_map_pos = np.array(gen_map_pos)
        gen_map_dist = np.roll(gen_map_pos, -1) - gen_map_pos
        gen_map_dist[-1] = np.nan
        #The minimum distance is 1e-10.
        low_dist_mask = gen_map_dist < 1e-10
        if np.sum(low_dist_mask) > 100:
            print "Warning: Over 100 genetic map distances were < 1e-10. These were changed to 1e-10. Probably unimportant, but noting as could e.g. indicate problem with genetic map."
        gen_map_dist[low_dist_mask] = 1e-10
    else:
        gen_map_dist = None
    print gen_map_dist
    
    curr_prob_low_archaic = copy.copy(init_prob_low_archaic)
    curr_prob_high_archaic = copy.copy(init_prob_high_archaic)
    viterbi_out_previous = np.zeros(len(string_of_archaic_alleles)) * np.nan
    for step in xrange(iterations):
        #Iteratively calculate the Viterbi algorithm, updating one or both of the prob_low_archaic and prob_high_archaic parameters.
        viterbi_out = support_viterbi(*support_construct_inputs(string_of_archaic_alleles, init_prob_high, curr_prob_low_archaic, curr_prob_high_archaic, transition_probability, gen_map_dist))
        suggested_curr_prob_low_archaic, suggested_curr_prob_high_archaic = support_re_calculate_freqs(viterbi_out, string_of_archaic_alleles, curr_prob_low_archaic, curr_prob_high_archaic)
        if fit_transitions is True and gen_map_dist is None:
            transition_probability = support_re_calculate_transition_probs(viterbi_out, string_of_archaic_alleles, transition_probability)
        elif fit_transitions is True and gen_map_dist is not None:
            transition_probability = support_re_calculate_transition_probs_genmap(viterbi_out, string_of_archaic_alleles, transition_probability, gen_map_dist)
        else:
            pass
        if 0 not in fixed_parameters:
            curr_prob_low_archaic = copy.copy(suggested_curr_prob_low_archaic)
        if 1 not in fixed_parameters:
            curr_prob_high_archaic = copy.copy(suggested_curr_prob_high_archaic)
        if np.sum(viterbi_out_previous == viterbi_out) == len(viterbi_out):
            print "Converged at step %d" %(step)
            break
        else:
            viterbi_out_previous = copy.copy(viterbi_out)
    return viterbi_out, string_of_archaic_alleles, [curr_prob_low_archaic, curr_prob_high_archaic], [np.sum(string_of_archaic_alleles == i) for i in range(max_str + 1)], transition_probability.flatten()

def support_viterbi(observation_space, state_space, initial_probabilities, sequence_of_observations, transition_matrix, emission_matrix, gen_map_dist):
    #observation_space of length N
    #state_space of length K
    #initial_probabilities of length K
    #sequence_of_observations of length T
    #transition_matrix of size K.K
    #emission_matrix of size K.N

    #The approach is to fill out two tables, starting with the initial state, given the probability of the observation in each state
    #For info see https://en.wikipedia.org/wiki/Viterbi_algorithm
    K = len(state_space) #num_states
    T = len(sequence_of_observations)
    z = np.ones(T, dtype = int) * -1
    X = np.ones(T, dtype = int) * -1
    
    T_1 = np.zeros((K, T))
    T_2 = np.zeros((K, T))
    for state in xrange(K):
        T_1[state, 0] = initial_probabilities[state] + emission_matrix[state,sequence_of_observations[0]]
    #print T_1
    #print T_2, "\n\n"
    for observation in xrange(1, T):
        
        for state in xrange(K):
            
            if gen_map_dist is None:
                #Constant transition probability per SNP, which is already log
                T_1[state, observation] = np.max(T_1[::, observation - 1] + transition_matrix[::,state]) + emission_matrix[state, sequence_of_observations[observation]]
                T_2[state, observation] = np.argmax(T_1[::, observation - 1] + transition_matrix[::,state])
            elif True:
                #Variable transition probability, depending on the genetic map. In this case, the transition matrix saves the switch rate per 1cM for 0->1 and 1->0. The converse is 1-this.
                #Transition rate is transitions per cM
                transition_probs  = np.array([0.0,0.0])
                alt_state = abs(state - 1)
                #RULE: at long distances, a shift may have a probability > 1 in the linear approximation. As this is not possible I set a maximum transition probability of 0.5.
                transition_probs[alt_state] = min([transition_matrix[::,state][alt_state] + np.log(gen_map_dist[observation - 1]), np.log(0.5)])
                transition_probs[state] = np.log(1.0 - np.exp(transition_probs[alt_state]))
                if np.sum(np.isnan(transition_probs)) > 0:
                    #Error, print some info
                    print T_1[::,observation - 1], T_1[::,observation]
                    print transition_matrix[::,state][alt_state], gen_map_dist[observation - 1], transition_probs, np.exp(transition_probs)
                    return None
                T_1[state, observation] = np.max(T_1[::, observation - 1] + transition_probs) + emission_matrix[state, sequence_of_observations[observation]]
                T_2[state, observation] = np.argmax(T_1[::, observation - 1] + transition_probs)
            else:
                #Placeholder. A more accurate rate would be p = m(1-exp[-rt]) ? Not tried to implement yet.
                pass
    #Go backwards down the possible paths (based on the most probable final state)
    z[T - 1] = np.argmax(T_1[::, T - 1]) #which is the more probable final state
    X[T - 1] = state_space[z[T - 1]]
    for observation in xrange(T - 1, 0, -1):
        z[observation - 1] = T_2[z[observation], observation] #The state at t is set to T_2[state at t+1, t]
        X[observation - 1] = state_space[z[observation - 1]]
    return X

def support_construct_inputs(string_of_archaic_alleles, init_prob_high, init_prob_low_archaic, init_prob_high_archaic, transition_probability, gen_map_dist):
    observation_space = np.array([0,1])
    state_space = np.array([0,1])
    initial_probabilities = np.array([np.log(1.0 - np.exp(init_prob_high)), init_prob_high]) #probability to start in each state, given the population
    sequence_of_observations = string_of_archaic_alleles #just the sequences of *observed* 0/1s
    transition_matrix = np.array(transition_probability)
    emission_matrix = np.array([[np.log(1.0 - np.exp(init_prob_low_archaic)), init_prob_low_archaic],
                                [np.log(1.0 - np.exp(init_prob_high_archaic)), init_prob_high_archaic]]) #this is the probability of emissions in the two hidden states
    return observation_space, state_space, initial_probabilities, sequence_of_observations, transition_matrix, emission_matrix, gen_map_dist

def support_re_calculate_freqs(state_path, sequence_of_observations, curr_low_archaic, curr_high_archaic):
    #Calculate the freqs for high and low from the sequence of observations
    print "Current emission probs:", np.exp(curr_low_archaic), np.exp(curr_high_archaic)
    low_obs = sequence_of_observations[state_path == 0]
    high_obs = sequence_of_observations[state_path == 1]
    tot_archaic_low_obs = np.sum(low_obs)
    tot_archaic_high_obs = np.sum(high_obs)
    if len(low_obs) != 0 and tot_archaic_low_obs != 0.0:
        prob_low_archaic = np.log(tot_archaic_low_obs / float(len(low_obs)))
    else:
        prob_low_archaic = curr_low_archaic
    if len(high_obs) != 0 and tot_archaic_high_obs != 0.0:
        prob_high_archaic = np.log(tot_archaic_high_obs / float(len(high_obs)))
    else:
        prob_high_archaic = curr_high_archaic
    return prob_low_archaic, prob_high_archaic

def support_re_calculate_transition_probs(viterbi_out, string_of_archaic_alleles, transition_probability):
    #For each hidden state, the number of possible transitions is how many times it occured, -1 if it is the last state.
    #Count what actually happened.
    print "Current transition probs:", np.exp(transition_probability)
    suggested_transition_probability = copy.copy(transition_probability)
    for hidden_state in range(len(transition_probability)):
        observed_following_states = viterbi_out[1:][viterbi_out[0:-1] == hidden_state]
        num_trans_opportunities = len(observed_following_states)
        if num_trans_opportunities > 1:
            for to_state in range(len(transition_probability)):
                suggested_transition_probability[hidden_state][to_state] = np.sum(observed_following_states == to_state) / float(num_trans_opportunities)
        else:
            #No information, so don't update
            pass
    return np.log(suggested_transition_probability)

def support_re_calculate_transition_probs_genmap(viterbi_out, string_of_archaic_alleles, transition_probability, gen_map_dist):
    #For each hidden state, the number of possible transitions is how many times it occured, -1 if it is the last state.
    #Count what actually happened.
    #I am trying to fit 'recombinations per cM' of the two transitions [0->1] and [1->0] only. The diagonals don't matter.
    print "Current transition probs:", np.exp(transition_probability)
    suggested_transition_probability = copy.copy(transition_probability)
    for hidden_state in range(len(transition_probability)):
        viterbi_hidden = [viterbi_out[0:-1] == hidden_state]
        observed_following_states = viterbi_out[1:][viterbi_hidden]
        following_gmap = gen_map_dist[0:-1][viterbi_hidden] #NB observed transition is IDx + 1 while the genamp distance is IDx
        num_trans_opportunities = len(observed_following_states)
        total_recom_distance = np.sum(following_gmap)
        #print "Hidden State:", hidden_state, "\tNum transition opportunities", num_trans_opportunities, "\tTotal recombination distance", total_recom_distance
        if num_trans_opportunities > 1:
            alt_state = abs(hidden_state - 1)
            to_state_mask = observed_following_states == alt_state
            transitions_per_cm = np.sum(to_state_mask) / float(total_recom_distance)
            suggested_transition_probability[hidden_state][alt_state] = transitions_per_cm
        else:
            #No information, so don't update
            pass
        suggested_transition_probability[hidden_state][hidden_state] = 1.0
    return np.log(suggested_transition_probability)
            

def support_calculate_viterbi_probabililty(state_path, sequence_of_observations, prob_low_archaic, prob_high_archaic, transition_probability):
    #Added 08/11/2017 to calculate the overall probability of the obsservations given the inferred parameters and state path.
    #This adds - the number of transitions, the total amount of 0 and 1 in the high introgression class and the total amount of 0 and 1 in the low introgression class
    #A bigger negative indicates a lower probability of this combination of states and observations, given the parameters of the model.
    high_introgression_class = state_path == 1
    low_introgression_class = state_path == 0
    total_high_introgression_class = np.sum(high_introgression_class)
    total_low_introgression_class = np.sum(low_introgression_class)
    total_1_in_high_introgression_class = np.sum(sequence_of_observations[high_introgression_class] == 1)
    total_1_in_low_introgression_class = np.sum(sequence_of_observations[low_introgression_class] == 1)
    total_0_in_high_introgression_class = total_high_introgression_class - total_1_in_high_introgression_class
    total_0_in_low_introgression_class = total_low_introgression_class - total_1_in_low_introgression_class
    total_transitions = np.sum(state_path - np.roll(state_path, 1)) - (1 if state_path[0] != state_path[-1] else 0)
    following_low_prob = transition_probability[0][state_path[1:][low_introgression_class[0:-1]]]
    following_high_prob = transition_probability[1][state_path[1:][high_introgression_class[0:-1]]]
    total_prob = sum([total_0_in_low_introgression_class * np.log(1.0 - np.exp(prob_low_archaic)),
                      total_1_in_low_introgression_class * prob_low_archaic,
                      total_0_in_high_introgression_class * np.log(1.0 - np.exp(prob_high_archaic)),
                      total_0_in_high_introgression_class * prob_high_archaic,
                      np.sum(following_low_prob),
                      np.sum(following_high_prob)])
    return total_prob
