#!/usr/bin/env python2.7

# ---- MODULE DOCSTRING
from __future__ import print_function,division
__doc__ = """

(C) Hive, Romain Wuilbercq, 2017
     _
    /_/_      .'''.
 =O(_)))) ...'     `.
    \_\              `.    .'''X
                       `..'
.---.  .---..-./`) ,---.  ,---.   .-''-.
|   |  |_ _|\ .-.')|   /  |   | .'_ _   \
|   |  ( ' )/ `-' \|  |   |  .'/ ( ` )   '
|   '-(_{;}_)`-'`"`|  | _ |  |. (_ o _)  |
|      (_,_) .---. |  _( )_  ||  (_,_)___|
| _ _--.   | |   | \ (_ o._) /'  \   .---.
|( ' ) |   | |   |  \ (_,_) /  \  `-'    /
(_{;}_)|   | |   |   \     /    \       /
'(_,_) '---' '---'    `---`      `'-..-'

The Artificial Bee Colony (ABC) algorithm is based on the
intelligent foraging behaviour of honey bee swarm, and was first proposed
by Karaboga in 2005.

Author:
------

Romain Wuilbercq

"""

# ---- IMPORT MODULES
try:
	import numpy as np
except:
	raise ImportError("Numpy module not installed.")

import time
import math	
import random
import sys
import copy

globSequences = []
###################################
#Codigo de fatorial retirado de https://stackoverflow.com/questions/16325988/factorial-of-a-large-number-in-python
def range_prod(lo,hi):
    if lo+1 < hi:
        mid = (hi+lo)//2
        return range_prod(lo,mid) * range_prod(mid+1,hi)
    if lo == hi:
        return lo
    return lo*hi

def treefactorial(n):
    if n < 2:
        return 1
    return range_prod(1,n)
 ####################################   
def thresholdConsensus(dna_sequences, consensus):
	dna_approvedSequences = []
	seqId = 0
	for sequence in dna_sequences:
		similarityRate = 0 #contador de similaridade com o consenso para cada sequencia
		i = 0
		while i < len(consensus):
			if sequence[i] == consensus[i]:
				similarityRate +=1
			i += 1
		similarityRate = similarityRate/len(consensus) #average
		#print("Seq",seqId,"threshold approval% =",similarityRate)
		if similarityRate > 0.50: #aprovacao
			dna_approvedSequences.append(sequence)
		seqId += 1
	
	#print("#Approved :",len(dna_approvedSequences))
	return dna_approvedSequences
		
def consensusMotif(positionCountMatrix):
	consensus = []
	i = 0

	while i < len(positionCountMatrix[0]):
		w_consensus = [] #para evitar bias para determinada base
		window = [positionCountMatrix[0][i],positionCountMatrix[1][i],positionCountMatrix[2][i],positionCountMatrix[3][i]]
		if window[0] >= window[1] and window[0] >= window[2] and window[0] >= window[3]:
			w_consensus.append('A')
		if window[1] >= window[0] and window[1] >= window[2] and window[1] >= window[3]:
			w_consensus.append('C')
		if window[2] >= window[0] and window[2] >= window[1] and window[2] >= window[3]:
			w_consensus.append('G')
		if window[3] >= window[0] and window[3] >= window[1] and window[3] >= window[2]:
			w_consensus.append('T')
		
		consensus.append(random.choice(w_consensus)) # um dos empatantes eh escolhido com probabilidades equivalentes
		
		#print (consensus[i], end = '')
		i += 1
	#print('\n')		
	return consensus
"""
def printDNAMatrix(DNAMatrix):
	names = ['A','C','G','T']
	i = 0
	for baseVector in DNAMatrix:
		print("#",names[i],":",end=" ")
		for base in baseVector:
			print(format(round(base,3),'.2f'), end=" ")
		print('\n')
		i += 1
	return
"""

def positionCountMatrix(dna_subsequences): 
    #col_1 col_2 ... col_n
    #Linha A
    #Linha C
    #Linha G
    #Linha T
    #o tamanho de todas as sequencias eh igual
    sequenceSize = len(dna_subsequences[0])
    #print(dna_subsequences)
    positionCountMatrix = np.zeros([4,sequenceSize],dtype=int)
    #print(len(dna_subsequences))
    for j in range(len(dna_subsequences)):
        sequence = dna_subsequences[j]
        for i in range(len(sequence)):
            if sequence[i] == 'A':
                positionCountMatrix[0][i] += 1
            elif sequence[i] == 'C':
                positionCountMatrix[1][i] += 1
            elif sequence[i] == 'G':
                positionCountMatrix[2][i] += 1
            elif sequence[i] == 'T':
                positionCountMatrix[3][i] += 1
        
    #print ('PCM')
    #printDNAMatrix(positionCountMatrix)
    return positionCountMatrix

def evaluator(vector):
	"""
	A n-dimensional Rastrigin's function is defined as:

							n
			f(x) = 10*n + Sigma { x_i^2 - 10*cos(2*PI*x_i) }
						   i=1

	where  -5.12 <= x_i <= 5.12.

	Thus the global minima of the function being f(x) = 0 at all x_i = 0.

	"""

	vector = np.array(vector)

	return 10 * vector.size + sum(vector*vector - 10 * np.cos(2 * np.pi * vector))

def sequenceFromSolution(dna_sequences, solutionVector):
    dna_solutionInstance = []
    dna_subSequences = []
    sequenceSize = len(dna_sequences[0])
    lowerLimit = 7
    higherLimit = 64
    if sequenceSize < higherLimit: #ajusta o tamanho da sequencia para um menor valor
        higherLimit = sequenceSize
    if sequenceSize >= lowerLimit: #uma sequencia nao eh considerada como motivo se for menor que 7
        motifSize = solutionVector[0] #tamanho do motivo
        i = 1
        #motifSize = 7 #tamanho do motivo
        #print("candidate size = ",motifSize)
        for sequence in dna_sequences:
            motifStart = solutionVector[i] #-1 pois o vetor comeca no 0
            #print("start = ",motifStart)
            subSequence = sequence[motifStart:motifStart+motifSize]
            dna_subSequences.append(subSequence)
            #print(subSequence)
            i += 1
    else:
        print("Sequencias de tamanho insuficiente:",sequenceSize,"<",lowerLimit)
    return dna_subSequences
def readFasta(filePath):
	
	dna_sequences = []
	fasta_sequences = SeqIO.parse(open(filePath),'fasta')
	i = 0
	for fasta in fasta_sequences:
		name,sequence = fasta.id,str(fasta.seq)
		dna_sequences.append(sequence)
		#print(">Seq",i,dna_sequences[len(dna_sequences)-1])
		i += 1
	return dna_sequences
	
def randomSubSequences(dna_sequences):
    dna_sequences = dna_sequences
    dna_solutionInstance = []
    dna_subSequences = []
    sequenceSize = len(dna_sequences[0])
    lowerLimit = 7
    higherLimit = 64
    if sequenceSize < higherLimit: #ajusta o tamanho da sequencia para um menor valor
        higherLimit = sequenceSize

    if sequenceSize >= lowerLimit: #uma sequencia nao eh considerada como motivo se for menor que 7
        motifSize = random.randint(lowerLimit,higherLimit) #tamanho do motivo
        #motifSize = 7 #tamanho do motivo
        #print("candidate size = ",motifSize)
        dna_solutionInstance.append(motifSize)

        for sequence in dna_sequences:
            motifStart = random.randint(0,sequenceSize-motifSize) #-1 pois o vetor comeca no 0
            #print("start = ",motifStart)
            subSequence = sequence[motifStart:motifStart+motifSize]
            dna_subSequences.append(subSequence)
            #print(subSequence)
            dna_solutionInstance.append(motifStart)
        #print("Sequencias de tamanho insuficiente:",sequenceSize,"<",lowerLimit)

    returnValues = []
    returnValues.append(dna_solutionInstance)
    returnValues.append(dna_subSequences)
    return returnValues

def printDNAMatrix(DNAMatrix):
    names = ['A','C','G','T']
    i = 0
    for baseVector in DNAMatrix:
        print("#",names[i],":",end=" ")
        for base in baseVector:
            print(base, end=" ")
            #print(format(round(base,3),'.2f'), end=" ")
        print('\n')
        i += 1
    return
def positionFrequencyMatrix(positionCountMatrix): #calcula as frequencias
    i = 0
    sumWindow = 0.0
    window = []
    sequenceSize = len(positionCountMatrix[0])
    positionFrequencyMatrix = np.zeros([4,sequenceSize]) #cria uma matriz de floats
    for i in range(sequenceSize):
        window = [positionCountMatrix[0][i],positionCountMatrix[1][i],positionCountMatrix[2][i],positionCountMatrix[3][i]]
        sumWindow = window[0] + window[1] + window[2] + window[3] #somatorio de nucleotideos na posicao i
        positionFrequencyMatrix[0][i] = float(window[0]/sumWindow)
        positionFrequencyMatrix[1][i] = float(window[1]/sumWindow)
        positionFrequencyMatrix[2][i] = float(window[2]/sumWindow)
        positionFrequencyMatrix[3][i] = float(window[3]/sumWindow)
    #print ('PFM')
    #printDNAMatrix(positionFrequencyMatrix)
    return positionFrequencyMatrix

def similarity(positionFrequencyMatrix):
	i = 0   
	sequenceSize = len(positionFrequencyMatrix[0])
	maxSum = 0.0
	while i < sequenceSize:
		window = [positionFrequencyMatrix[0][i],positionFrequencyMatrix[1][i],positionFrequencyMatrix[2][i],positionFrequencyMatrix[3][i]]
		maxSum += max(window)
		i += 1

	similarity = maxSum/sequenceSize
	return similarity

def complexity(motif):
	motifSize = len(motif)

	numA = 0
	numC = 0
	numG = 0
	numT = 0
	motifSizeFactorial = treefactorial(motifSize)
	i = 0
	while i < motifSize: #conta a quantidade de cada base
		if motif[i] == 'A':
			numA += 1
				
		elif motif[i] == 'C':
			numC += 1

		elif motif[i] == 'G':
			numG += 1

		elif motif[i] == 'T':
			numT += 1
		i += 1

	productBasesFactorial = treefactorial(numA)*treefactorial(numC)*treefactorial(numG)*treefactorial(numT)
							#produtorio do fatorial do numero de cada base

	complexity = math.log(motifSizeFactorial/productBasesFactorial,4) #logaritmo base 4 (numero de bases) 
	return complexity

def isBiased(support,totalSequences):
	if totalSequences <= 4 :#minimo dado pelos autores (Alvarez)
		if support >= 2 :
			isBiased = False
		else :
			isBiased = True
	else :
		if support >=3 :
			isBiased = False
		else :
			isBiased = True

	return isBiased
	

	



# ---- BEE CLASS
class CandidateMotif(object):
    # representacao do motivo candidato
    def __init__(self,motif,solution,similarity,complexity,support,consensus):
        self.motif = motif
        self.solution = solution
        self.similarity = similarity
        self.support = support
        self.complexity = complexity
        self.consensus = consensus

    def _printSolution(self,sequence):
        finalSubSequences = sequenceFromSolution(sequence,self.solution)         
        dna_approvedSequences = thresholdConsensus(finalSubSequences,self.consensus)
        
        negSolution = self._negativateSolution(sequence)
        
        
        print("=========================")
        print("Final motif:",end = "")
        for base in self.motif:
            print(base, end = "")
        print("")
        """
        print(self.solution)
        print(negSolution)
        """
        i = 1
        while i < len(self.solution):
            currSeq = sequence[i-1]
            substring = currSeq[self.solution[i]:self.solution[i]+self.solution[0]]
            for approved in dna_approvedSequences:  
                #print(approved,"?")
                if approved == substring:
                    #print(approved)
                    print(i-1,",",negSolution[i],",",substring, sep='') #[self.solution[i]:self.solution[i]+self.solution[0]]
                    #print(i-1,",",negSolution[i],",",currSeq[negSolution[i]:negSolution[i]+self.solution[0]-1], sep='') #[self.solution[i]:self.solution[i]+self.solution[0]]
                    #print(negSolution[i],":",negSolution[i]+self.solution[0]-1)
            i += 1
        print("Support", self.support)
        print("Similarity",self.similarity)
        print("Complexity",self.complexity)
        print("=========================")
    
    def _negativateSolution(self,sequence):
        from copy import deepcopy
        negSolution = deepcopy(self.solution)
        seqLen = len(sequence[0])
        i = 1
        while i < len(self.solution):
            negSolution[i]-=seqLen
            i += 1
        
        return negSolution 

class Bee(object):
    """ Creates a bee object. """

    def __init__(self, lower, upper, dna_sequences,obj, funcon=None):
        """

        Instantiates a bee object randomly.

        Parameters:
        ----------
            :param list lower  : lower bound of solution vector
            :param list upper  : upper bound of solution vector
            :param def  fun    : evaluation function
            :param def  funcon : constraints function, must return a boolean

        """
        self.obj = obj
        self.dna_sequences = dna_sequences
        #print(dna_sequences)
        biased = True        #enquanto nao houver solucao valida, instancia uma nova
        resultVector = []
        solution = -1
        self.dna_subSequence = []
        consensus = []
        pcm = []
        motifSupport = []
        dna_approvedSequences = []

        while biased == True:
            resultVector = randomSubSequences(self.dna_sequences)
            solution = resultVector[0]
            self.dna_subSequences = resultVector[1]
            pcm = positionCountMatrix(self.dna_subSequences)
            #print(self.dna_subSequences)
            #printDNAMatrix(pcm)
            consensus = consensusMotif(pcm)
            
            dna_approvedSequences = thresholdConsensus(self.dna_subSequences,consensus)
            
            motifSupport= len(dna_approvedSequences)
            biased = isBiased(motifSupport,len(self.dna_subSequences))
        
        finalPcm = positionCountMatrix(dna_approvedSequences)
        finalPfm = positionFrequencyMatrix(finalPcm)
        finalMotif = consensusMotif(finalPcm)
        motifSimilarity = similarity(finalPfm)
        motifComplexity = complexity(finalMotif)
        
        self.solutionVector = solution	    
        self.candidate = CandidateMotif(finalMotif,solution,motifSimilarity,motifComplexity,motifSupport,consensus)
        self.vector = solution
        self.valid = biased
        #self.candidate._printSolution()

        # creates a random solution vector
        #self._random(lower, upper)
        """
        # checks if the problem constraint(s) are satisfied
        if not funcon:
            self.valid = True
        else:
            self.valid = funcon(self.vector)

        # computes fitness of solution vector
        if (fun != None):
            self.value = fun(self.vector)
        else:
            self.value = sys.float_info.max
        
        """
        if(self.obj == 'similarity'):
            self.value = self.candidate.similarity
        elif(self.obj == 'support'):
            self.value = self.candidate.support
        elif(self.obj == 'complexity'):
            self.value = self.candidate.complexity
        
        self._fitness()
        # initialises trial limit counter - i.e. abandonment counter
        
        self.counter = 0
    

        


    def _random(self, lower, upper):
        """ Initialises a solution vector randomly. """

        self.vector = []
        for i in range(len(lower)):
            self.vector.append( lower[i] + random.random() * (upper[i] - lower[i]) )

    def _fitness(self):
        """

        Evaluates the fitness of a solution vector.

        The fitness is a measure of the quality of a solution.

        """
        """
        if (self.value >= 0):
            self.fitness = 1 / (1 + self.value)
        else:
            self.fitness = 1 + abs(self.value)
        """
        self.fitness = self.value
class BeeHive(object):
    """

    Creates an Artificial Bee Colony (ABC) algorithm.

    The population of the hive is composed of three distinct types
    of individuals:

        1. "employees",
        2. "onlookers",
        3. "scouts".

    The employed bees and onlooker bees exploit the nectar
    sources around the hive - i.e. exploitation phase - while the
    scouts explore the solution domain - i.e. exploration phase.

    The number of nectar sources around the hive is equal to
    the number of actively employed bees and the number of employees
    is equal to the number of onlooker bees.

    """

    def run(self):
        """ Runs an Artificial Bee Colony (ABC) algorithm. """

        cost = {}; cost["best"] = []; cost["mean"] = []
        
        for itr in range(self.max_itrs):

            # employees phase
            for index in range(self.size):
                self.send_employee(index)
            
            # onlookers phase
            self.send_onlookers()

            # scouts phase
            self.send_scout()

            # computes best path
            self.find_best()

            # stores convergence information
            cost["best"].append( self.best )
            cost["mean"].append( sum( [ bee.value for bee in self.population ] ) / self.size )

            # prints out information about computation
            if self.verbose:
                self._verbose(itr, cost)
        #print("\n")
        #print("****BEST SOLUTION****")
        #self.bestSolution._printSolution()

        #sequenceFromSolution(self.dna_sequences,self.bestSolution.solution)
            
        #for best in cost["best"]:
            #print(best)
        return cost

    def __init__(self                 ,
                 lower, upper, dna_sequences,
                 numb_bees, max_itrs, max_trials,obj,
                 fun          = None  ,
                 selfun       = None  ,
                 seed         = None  ,
                 verbose      = False ,
                 extra_params = None  ,):
        """

        Instantiates a bee hive object.

        1. INITIALISATION PHASE.
        -----------------------

        The initial population of bees should cover the entire search space as
        much as possible by randomizing individuals within the search
        space constrained by the prescribed lower and upper bounds.

        Parameters:
        ----------

            :param list lower          : lower bound of solution vector
            :param list upper          : upper bound of solution vector
            :param def fun             : evaluation function of the optimal problem
            :param def numb_bees       : number of active bees within the hive
            :param int max_trials      : max number of trials without any improvment
            :param def selfun          : custom selection function
            :param int seed            : seed of random number generator
            :param boolean verbose     : makes computation verbose
            :param dict extra_params   : optional extra arguments for selection function selfun

        """

        # checks input
        assert (len(upper) == len(lower)), "'lower' and 'upper' must be a list of the same length."

        """# generates a seed for the random number generator
        if (seed == None):
            self.seed = random.randint(0, 1000)
        else:
            self.seed = seed
        random.seed(self.seed)
        """
        self.obj = obj
        self.dna_sequences = dna_sequences
        #print("AAA",len(self.dna_sequences))
        
        #if(len(self.dna_sequences)<3):
            #
        # computes the number of employees
        self.size = int((numb_bees + numb_bees % 2))

        # assigns properties of algorithm
        self.dim = len(dna_sequences)+1 
        self.max_itrs = max_itrs
        #if (max_trials == None):
        #    self.max_trials = 0.6 * self.size * self.dim
        #else:
        #    self.max_trials = max_trials
        self.max_trials = max_trials
        self.selfun = selfun
        self.extra_params = extra_params

        # assigns properties of the optimisation problem
        self.evaluate = fun
        self.lower    = lower
        self.upper    = upper

        # initialises current best and its a solution vector
        self.best = 0
        self.solution = None

        # creates a bee hive
      
        #sequenceSize-motifSize !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        self.population = [ Bee(lower, upper, self.dna_sequences,self.obj) for i in range(self.size) ]

        # initialises best solution vector to food nectar
        self.find_best()

        # computes selection probability
        self.compute_probability()

        # verbosity of computation
        self.verbose = verbose

    def find_best(self):
        """ Finds current best bee candidate. """

        values = [ bee.value for bee in self.population ]
        index  = values.index(max(values)) #Maximize
        if (values[index] > self.best):
            self.best     = values[index]
            self.bestSolution = self.population[index].candidate

    def compute_probability(self):
        """

        Computes the relative chance that a given solution vector is
        chosen by an onlooker bee after the Waggle dance ceremony when
        employed bees are back within the hive.

        """

        # retrieves fitness of bees within the hive
        values = [bee.fitness for bee in self.population]
        max_values = max(values)

        # computes probalities the way Karaboga does in his classic ABC implementation
        if (self.selfun == None):
            self.probas = [0.9 * v / max_values + 0.1 for v in values]
        else:
            if (self.extra_params != None):
                self.probas = self.selfun(list(values), **self.extra_params)
            else:
                self.probas = self.selfun(values)

        # returns intervals of probabilities
        return [sum(self.probas[:i+1]) for i in range(self.size)]

    def send_employee(self, index):
        """

        2. SEND EMPLOYED BEES PHASE.
        ---------------------------

        During this 2nd phase, new candidate solutions are produced for
        each employed bee by cross-over and mutation of the employees.

        If the modified vector of the mutant bee solution is better than
        that of the original bee, the new vector is assigned to the bee.

        """

        # deepcopies current bee solution vector
        zombee = copy.deepcopy(self.population[index])

        # draws a dimension to be crossed-over and mutated
        d = random.randint(0, self.dim-1)

        # selects another bee
        bee_ix = index;
        bee_ix = index;
        while (bee_ix == index): bee_ix = random.randint(0, self.size-1)

        # produces a mutant based on current bee and bee's friend
        zombee.vector[d] = self._mutate(d, index, bee_ix)
        
        # computes fitness of mutant
        if(zombee.obj == 'similarity'):
            zombee.value = zombee.candidate.similarity
        elif(zombee.obj == 'support'):
            zombee.value = zombee.candidate.support
        elif(zombee.obj == 'complexity'):
            zombee.value = zombee.candidate.complexity
        zombee._fitness()

        # deterministic crowding
        if (zombee.fitness > self.population[index].fitness):
            self.population[index] = copy.deepcopy(zombee)
            self.population[index].counter = 0
        else:
            self.population[index].counter += 1
        

    def send_onlookers(self):
        """

        3. SEND ONLOOKERS PHASE.
        -----------------------

        We define as many onlooker bees as there are employed bees in
        the hive since onlooker bees will attempt to locally improve the
        solution path of the employed bee they have decided to follow
        after the waggle dance phase.

        If they improve it, they will communicate their findings to the bee
        they initially watched "waggle dancing".

        """

        # sends onlookers
        numb_onlookers = 0; beta = 0
        while (numb_onlookers < self.size):

            # draws a random number from U[0,1]
            phi = random.random()

            # increments roulette wheel parameter beta
            beta += phi * max(self.probas)
            beta %= max(self.probas)

            # selects a new onlooker based on waggle dance
            index = self.select(beta)

            # sends new onlooker
            self.send_employee(index)

            # increments number of onlookers
            numb_onlookers += 1

    def select(self, beta):
        """

        4. WAGGLE DANCE PHASE.
        ---------------------

        During this 4th phase, onlooker bees are recruited using a roulette
        wheel selection.

        This phase represents the "waggle dance" of honey bees (i.e. figure-
        eight dance). By performing this dance, successful foragers
        (i.e. "employed" bees) can share, with other members of the
        colony, information about the direction and distance to patches of
        flowers yielding nectar and pollen, to water sources, or to new
        nest-site locations.

        During the recruitment, the bee colony is re-sampled in order to mostly
        keep, within the hive, the solution vector of employed bees that have a
        good fitness as well as a small number of bees with lower fitnesses to
        enforce diversity.

        Parameter(s):
        ------------
            :param float beta : "roulette wheel selection" parameter - i.e. 0 <= beta <= max(probas)

        """

        # computes probability intervals "online" - i.e. re-computed after each onlooker
        probas = self.compute_probability()

        # selects a new potential "onlooker" bee
        for index in range(self.size):
            if (beta < probas[index]):
                return index

    def send_scout(self):
        """

        5. SEND SCOUT BEE PHASE.
        -----------------------

        Identifies bees whose abandonment counts exceed preset trials limit,
        abandons it and creates a new random bee to explore new random area
        of the domain space.

        In real life, after the depletion of a food nectar source, a bee moves
        on to other food sources.

        By this means, the employed bee which cannot improve their solution
        until the abandonment counter reaches the limit of trials becomes a
        scout bee. Therefore, scout bees in ABC algorithm prevent stagnation
        of employed bee population.

        Intuitively, this method provides an easy means to overcome any local
        optima within which a bee may have been trapped.

        """

        # retrieves the number of trials for all bees
        trials = [ self.population[i].counter for i in range(self.size) ]

        # identifies the bee with the greatest number of trials
        index = trials.index(max(trials))

        # checks if its number of trials exceeds the pre-set maximum number of trials
        if (trials[index] > self.max_trials):

            # creates a new scout bee randomly
            self.population[index] = Bee(self.lower, self.upper, self.dna_sequences,self.obj)

            # sends scout bee to exploit its solution vector
            self.send_employee(index)

    def _mutate(self, dim, current_bee, other_bee):
        """

        Mutates a given solution vector - i.e. for continuous
        real-values.

        Parameters:
        ----------

            :param int dim         : vector's dimension to be mutated
            :param int current_bee : index of current bee
            :param int other_bee   : index of another bee to cross-over

        """
        return self.population[current_bee].vector[dim]    + \
               random.choice([-1,1])                       * \
               (self.population[current_bee].vector[dim] - self.population[other_bee].vector[dim])


    def _verbose(self, itr, cost):
        """ Displays information about computation. """

        msg = "# Iter = {} | Best Evaluation Value = {} | Mean Evaluation Value = {} "
        print(msg.format(int(itr), cost["best"][itr], cost["mean"][itr]))

# ---- END
