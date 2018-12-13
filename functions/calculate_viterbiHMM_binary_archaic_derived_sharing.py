###This is an attempt to implement a toy HMM for detecting archaic admixture
###The method imagines re-implementing http://biorxiv.org/content/biorxiv/early/2015/07/15/022632.full.pdf
###It is actually apparently successful!
###The method will work to categorise the genome into chunks of two types, given:
###A binary state vector
###An unknown probability associated with each state for each type
###An fixed penalty for transitions between types
###No explicit penalty for e.g. longer chunks, or prior on the size of chunks,
###though the total number of transitions is impacted by the transition penalty
###
###There is potential for improvement when there is a strong prior on the chunk length.
###
###Implemented summer 2017.
###Modified 08/11/2017 to allow fixing of parameters. Also added code for calculating the probability of the inferred state.
###Modified 17/11/2017 to allow fitting of, and fixing of, the transition parameters too. The program also now expects the full transition matrix to be provided, allowing for unequal transition rates.
###Modified 19/11/2017 to allow the input of a genetic map

import numpy as np
import copy

"""
#string_of_archaic_alleles = np.array([0,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0])
string_of_archaic_alleles = []
for i in range(500): #500 is 500000 SNPs. We might actually have about 500k SNPs per chromosome, so this is probably about right...
    #Generate according to the model - windows of length 100 with probability
    if np.random.rand() < 0.05:
        string_of_archaic_alleles = string_of_archaic_alleles + [np.random.rand() < 0.5 for j in range(100)]
    else:
        string_of_archaic_alleles = string_of_archaic_alleles + [np.random.rand() < 0.05 for j in range(100)]
string_of_archaic_alleles = np.array(string_of_archaic_alleles, dtype = int)
#string should be 50k long, with 5% high density 100 SNP regions and 95% low density 100SNP regions
#string_of_archaic_alleles = np.array([0,1,1,0,0,1,0,0])
init_prob_high_PNG = np.log(0.5)
init_prob_high_Australian = np.log(0.5)
init_prob_high_archaic = np.log(0.0276)
init_prob_low_archaic = np.log(0.0058)
transition_probability = np.log(0.3) #-8.0
#initial_states = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1])
"""

def calculate_viterbiHMM_binary_archaic_derived_sharing(infile_types, init_prob_high, init_prob_low_archaic, init_prob_high_archaic, transition_probability, iterations = 100, fixed_parameters = [], fit_transitions = False, genetic_map = None):
    """
    This function calculates the viterbiHMM until it converges.
    Read in the input string file.
    Using the initial frequencies of high-admixture regions and the initial proportion of 1 alleles in low and high region, calculate the Viterbi algorithm.
    Use the output to re-calculate the frequencies of the different regions and the probabilities associated with them.
    Feed these probabilities back into the Viterbi algorithm and repeat until it stabilises.
    The process can be repeated multiple times with different parameters by inputting len>1 lists of values init_prob_high, init_prob_low_archaic, init_prob_high_archaic
    In this case, the iteration that has the highest log-probability is used.
    fixed_parameters specifies if any of the parameters are fixed at their initial value. In this case, 0 is for init_prob_low_archaic and 1 is for init_prob_high_archaic.
    There is an option to fit the transition probability using fit_transitions = True
    There is finally an option to incorporate a genetic map using genetic_map. The format should be one line per positions - position \t genetic map position
    """
    open_infile = gzip.open if infile_types[-3:] == '.gz' else open
    
    with open_infile(infile_types, 'rb') as f_in:
        #Note that 0-9 are characters and >9 are ascii to let me use up to 90-48 = 42 classes
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
        #The minimum distance is 1e-10. This is important.
        low_dist_mask = gen_map_dist < 1e-10
        if np.sum(low_dist_mask) > 100:
            print "Warning: Over 100 genetic map distances were < 1e-10. These were changed to 1e-10. Noting as could indicate problem with genetic map."
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
    #Table T_1 stores the probability of the most likely path so far
    #Table T_2 stores ???the most probable previous state given the current state???
    #I don't *quite* get it but this seems to work fine vs a brute force method explicitly recording each possible path!
    #Algorithm is from Wikipedia.
    K = len(state_space) #num_states
    T = len(sequence_of_observations)
    z = np.ones(T, dtype = int) * -1
    X = np.ones(T, dtype = int) * -1
    
    T_1 = np.zeros((K, T)) #stores the probability of the most likely path so far
    T_2 = np.zeros((K, T)) #stores the most probable previous state, given the current state????
    for state in xrange(K):
        T_1[state, 0] = initial_probabilities[state] + emission_matrix[state,sequence_of_observations[0]]
        #T_2[state, 0] = 0.0 #not needed, I think this isn't relevant as it would be tell us the t = -1 state
    #print T_1
    #print T_2, "\n\n"
    for observation in xrange(1, T):
        
        for state in xrange(K):
            
            #Set the next observation in T_1 to i) the probability of being
            #print state, observation, T_1[::, observation - 1] + transition_matrix[state], emission_matrix[state, sequence_of_observations[observation]]
            #Next T_1 is set to the log-probability of the observation for the suggested state plus the more probable of (previous state is current state and no transition) and (previous state isn't current state and transition)
            #And set T_2 to 0 if (previous_state is current state and no transition) is more likely
            #Or 1 if (previous_state is not current state and transition) is more likely
            #'Given that my current state is state, is it more likely that my previous state was 0 or 1?'
            if gen_map_dist is None:
                #Constant transition probability, which is already in log format
                T_1[state, observation] = np.max(T_1[::, observation - 1] + transition_matrix[::,state]) + emission_matrix[state, sequence_of_observations[observation]]
                T_2[state, observation] = np.argmax(T_1[::, observation - 1] + transition_matrix[::,state])
            elif True:
                #Variable transition probability, depending on the genetic map. In this case, the transition matrix save the switch rate per 1cM for 0->1 and 1->0. The converse is 1-this.
                #The linear equation, p(0->1) = rm(t-1) says that the transition matrix implicitly gives m(t-1) and (1-m)(t-1); but this will only be true when the tracts are being accurately estimated.
                #I prefer not to rely on this, but it is an optional constraint in the fitting.
                #This will be more tricky with multiple states.
                transition_probs  = np.array([0.0,0.0])
                alt_state = abs(state - 1)
                #RULE: at long distances, a shift may have a probability > 1 in the linear approximation.
                #Realistically, the shift most likely has probability m or (1-m). 0.5 is saying that the algorithm doesn't care what the previous state was.
                transition_probs[alt_state] = min([transition_matrix[::,state][alt_state] + np.log(gen_map_dist[observation - 1]), np.log(0.5)])
                transition_probs[state] = np.log(1.0 - np.exp(transition_probs[alt_state]))
                #print transition_matrix[::,state], np.shape(transition_matrix[::,state])
                if np.sum(np.isnan(transition_probs)) > 0:
                    print T_1[::,observation - 1], T_1[::,observation]
                    print transition_matrix[::,state][alt_state], gen_map_dist[observation - 1], transition_probs, np.exp(transition_probs)
                    return None
                T_1[state, observation] = np.max(T_1[::, observation - 1] + transition_probs) + emission_matrix[state, sequence_of_observations[observation]]
                T_2[state, observation] = np.argmax(T_1[::, observation - 1] + transition_probs)
            else:
                #Placeholder. This is a more accurate rate, p = m(1-exp[-rt]). It's not linear, so more annoying to implement.
                pass
            #print T_1[::,observation - 1], T_1[::,observation]
        #print T_1
        #print T_2, "\n"
    #Go backwards down the possible paths (based on the most probable final state)
    z[T - 1] = np.argmax(T_1[::, T - 1]) #which is the more probable final state
    X[T - 1] = state_space[z[T - 1]]
    for observation in xrange(T - 1, 0, -1):
        #print "Setting", observation - 1
        #print "Previous obsevation", z[observation]
        #print observation, z[observation]
        #print X
        z[observation - 1] = T_2[z[observation], observation] #The state at t is set to T_2[state at t+1, t]
        X[observation - 1] = state_space[z[observation - 1]]
        #print X, "\n"
    return X#, T_1, T_2

def support_construct_inputs(string_of_archaic_alleles, init_prob_high, init_prob_low_archaic, init_prob_high_archaic, transition_probability, gen_map_dist):
    observation_space = np.array([0,1])
    state_space = np.array([0,1])
    initial_probabilities = np.array([np.log(1.0 - np.exp(init_prob_high)), init_prob_high]) #probability to start in each state, given the population
    sequence_of_observations = string_of_archaic_alleles #just the sequences of 0/1 giving (not-E-allele)/(E-allele)
    transition_matrix = np.array(transition_probability)
    #transition_matrix = np.array([[np.log(1.0 - np.exp(transition_probability)), transition_probability],
    #                              [transition_probability, np.log(1.0 - np.exp(transition_probability))]]) #this is the transition between states
    emission_matrix = np.array([[np.log(1.0 - np.exp(init_prob_low_archaic)), init_prob_low_archaic],
                                [np.log(1.0 - np.exp(init_prob_high_archaic)), init_prob_high_archaic]]) #this is the probability of E-alleles in each state
    return observation_space, state_space, initial_probabilities, sequence_of_observations, transition_matrix, emission_matrix, gen_map_dist

def support_re_calculate_freqs(state_path, sequence_of_observations, curr_low_archaic, curr_high_archaic):
    #Calculate the freqs for high and low from the sequence of observations
    #Do I work out the overall frequency of frequency per window? Probably overall?
    print "Current emission probs:", np.exp(curr_low_archaic), np.exp(curr_high_archaic)
    low_obs = sequence_of_observations[state_path == 0]
    high_obs = sequence_of_observations[state_path == 1]
    tot_archaic_low_obs = np.sum(low_obs)
    tot_archaic_high_obs = np.sum(high_obs)
    print tot_archaic_low_obs, tot_archaic_high_obs
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
    #There is an option to constrain the rates according to m(t-1) and (1-m)(t-1).
    #This is more complicated, and requires me to explore different parameters.
    print "Current transition probs:", np.exp(transition_probability)
    suggested_transition_probability = copy.copy(transition_probability)
    for hidden_state in range(len(transition_probability)):
        viterbi_hidden = [viterbi_out[0:-1] == hidden_state]
        observed_following_states = viterbi_out[1:][viterbi_hidden]
        following_gmap = gen_map_dist[0:-1][viterbi_hidden] #Note that the observed transition is IDx + 1 while the genamp distance is IDx
        num_trans_opportunities = len(observed_following_states)
        total_recom_distance = np.sum(following_gmap) #This is the important one.
        print hidden_state, num_trans_opportunities, total_recom_distance
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

"""
###How about I do it myself. At each iteration, I add up all the probabilities. I do it on a small run
def brute_prob(observation_space, state_space, initial_probabilities, sequence_of_observations, transition_matrix, emission_matrix):
    probs_so_far = {}
    probs_so_far[0] = {}
    probs_so_far[0]['0'] = initial_probabilities[0] + emission_matrix[0,sequence_of_observations[0]]
    probs_so_far[0]['1'] = initial_probabilities[1] + emission_matrix[1,sequence_of_observations[0]]
    #print probs_so_far
    for i in range(len(sequence_of_observations) - 1):
        probs_so_far[i + 1] = {}
        for key in probs_so_far[i]:
            for state in [0,1]:
                new_key = key + str(state)
                previous_state = key[-1]
                probs_so_far[i + 1][new_key] = probs_so_far[i][key] + transition_matrix[int(previous_state), state] + emission_matrix[state,sequence_of_observations[i + 1]]
    most_probable_run = None
    most_probable_prob = None
    for j in probs_so_far[len(sequence_of_observations) - 1].items():
        if most_probable_run is None:
            most_probable_run = j[0]
            most_probable_prob = j[1]
        elif j[1] > most_probable_prob:
            most_probable_run = j[0]
            most_probable_prob = j[1]
        else:
            pass
    most_probable_run = np.array([i for i in most_probable_run], dtype = float)
    return most_probable_run #probs_so_far[i + 1],

def test():
    for tp in np.linspace(0.01, 0.99, 10):
        print tp
        for h in np.linspace(0.01,0.99, 10):
            for l in np.linspace(0.001, h, 10):
                inputs = support_construct_inputs(string_of_archaic_alleles, np.log(0.7), np.log(l), np.log(h), np.log(tp))
                a = support_viterbi(*inputs)
                b = brute_prob(*inputs)
                #print a, b
                if np.sum(a == b) == 8:
                    pass
                else:
                    print "DISAGREEMENT", tp, h, l

"""
