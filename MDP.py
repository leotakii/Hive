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

"""

# ---- IMPORT MODULES

try:
	import numpy as np
except:
	raise ImportError("Numpy module not installed.")

import random
import time
import math	
from Bio import SeqIO
#from Hive import DNASequence as dnaSeq
from Hive import Hive
from Hive import Utilities

random.seed(time.time()) #testando deterministicamente

# ---- CREATE TEST CASE
class CandidateMotif:
	# representacao do motivo candidato
	def __init__(self,vector):
		self.motif = finalMotif(vector)
		self.similarity = similarity(vector)
		self.support = support(vector)
		self.complexity = complexity(complexity)





def finalMotif():
	finalMotif = []
	
	candidate = CandidateMotif();
	

	return finalMotif



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
		print("Seq",seqId,"similarity % =",similarityRate)
		if similarityRate > 0.50: #aprovacao
			dna_approvedSequences.append(sequence)
		seqId += 1
	
	print("#Approved :",len(dna_approvedSequences))
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
		
		print (consensus[i], end = '')
		i += 1
	print('\n')		
	return consensus

def printDNAMatrix(DNAMatrix):
	names = ['A','C','G','T']
	i = 0
	for baseVector in DNAMatrix:
		print("#",names[i],":",end=" ")
		for base in baseVector:
			print(base, end=" ")
		print('\n')
		i += 1
	return

def positionCountMatrix(dna_sequences):
	sequenceSize = len(dna_sequences[0])
	positionCountMatrix = np.zeros([4,sequenceSize],dtype=int) #o tamanho de todas as sequencias eh igual	
														#col_1 col_2 ... col_n
											  #Linha A
											  #Linha C
											  #Linha G
											  #Linha T

	for sequence in dna_sequences:
		i = 0
		while i < len(sequence): #itera bases da sequencia
			if sequence[i] == 'A':
				positionCountMatrix[0][i] += 1
				
			elif sequence[i] == 'C':
				positionCountMatrix[1][i] += 1

			elif sequence[i] == 'G':
				positionCountMatrix[2][i] += 1

			elif sequence[i] == 'T':
				positionCountMatrix[3][i] += 1
			i += 1
	print ('PCM')
	printDNAMatrix(positionCountMatrix)
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





def readFasta(filePath):
	
	dna_sequences = []
	fasta_sequences = SeqIO.parse(open(filePath),'fasta')

	for fasta in fasta_sequences:
		name,sequence = fasta.id,str(fasta.seq)
		dna_sequences.append(sequence)
		print(dna_sequences[len(dna_sequences)-1]+'\n')
	
	return dna_sequences
	
def subSequences(dna_sequences):
	dna_subSequences = []

	motifSize = random.randint(7,64) #tamanho do motivo
	#motifSize = 7 #tamanho do motivo
	print("size = ",motifSize)
	sequenceSize = len(dna_sequences[0])
	
	for sequence in dna_sequences:
		motifStart = random.randint(0,sequenceSize-motifSize)
		print("start = ",motifStart)
		subSequence = sequence[motifStart:motifStart+motifSize]
		dna_subSequences.append(subSequence)
		print(subSequence)
		
	return dna_subSequences

def positionFrequencyMatrix(positionCountMatrix): #calcula as frequencias
	i = 0
	sumWindow = 0.0
	window = []
	sequenceSize = len(positionCountMatrix[0])
	positionFrequencyMatrix = np.zeros([4,sequenceSize]) #cria uma matriz de floats
	
	while i < sequenceSize:

		window = [positionCountMatrix[0][i],positionCountMatrix[1][i],positionCountMatrix[2][i],positionCountMatrix[3][i]]
		sumWindow = window[0] + window[1] + window[2] + window[3] #somatorio de nucleotideos na posicao i
		positionFrequencyMatrix[0][i] = float(window[0]/sumWindow)
		positionFrequencyMatrix[1][i] = float(window[1]/sumWindow)
		positionFrequencyMatrix[2][i] = float(window[2]/sumWindow)
		positionFrequencyMatrix[3][i] = float(window[3]/sumWindow)

		i += 1

	print ('PFM')
	printDNAMatrix(positionFrequencyMatrix)
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

def run():
	
	# creates model
	#ndim = int(10)
	#model = Hive.BeeHive(lower = [-5.12]*ndim  ,
	#                     upper = [ 5.12]*ndim  ,
	#                     fun       = evaluator ,
	#                     numb_bees =  50       ,
	#                     max_itrs  =  100       ,)

	# runs model
	#cost = model.run()

	# plots convergence
	#Utilities.ConvergencePlot(cost)

	# prints out best solution
	#print("Fitness Value ABC: {0}".format(model.best))
	



	dna_sequences = readFasta("assets/Real/testDNA.fasta")
	dna_subSequences = subSequences(dna_sequences)
	""" #bloco com sequencia inteira
	pcm = positionCountMatrix(dna_sequences)
	consensus = consensusMotif(pcm)
	thresholdConsensus(dna_sequences,consensus)
	"""
	pcm = positionCountMatrix(dna_subSequences)
	consensus = consensusMotif(pcm)
	dna_approvedSequences = thresholdConsensus(dna_subSequences,consensus)
	finalPcm = positionCountMatrix(dna_approvedSequences)
	finalPfm = positionFrequencyMatrix(finalPcm)
	finalMotif = consensusMotif(finalPcm)
	similarityMotif = similarity(finalPfm)
	complexityMotif = complexity(finalMotif)

	supportMotif = len(dna_approvedSequences)
	print("Similarity:",similarityMotif)
	print("Complexity:",complexityMotif)
	print("Support:",supportMotif)
	###############readFasta(sys.argv[1]) #1 = path name




	"""dna_sequences = []
	fasta_sequences = SeqIO.parse(open("assets/Real/dm01r.fasta"),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		#new_sequence = some_function(sequence)
		dna_sequences.append(sequence)
		print(dna_sequences[len(dna_sequences)-1]+'\n')
	"""
	#print("aa");

if __name__ == "__main__":
	run()


# ---- END
