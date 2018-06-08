#!/usr/bin/env python2.7

# ---- MODULE DOCSTRING

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

from Bio import SeqIO
#from Hive import DNASequence as dnaSeq
from Hive import Hive
from Hive import Utilities

# ---- CREATE TEST CASE
class CandidateMotif:
	# representacao do motivo candidato
	def __init__(self,vector):
		self.motif = finalMotif(vector)
		self.similarity = similarity(vector)
		self.support = support(vector)
		self.complexity = complexity(complexity)






def positionCountMatrix(dna_sequences):
	positionCountMatrix = np.zeros([4,len(dna_sequences[0])],dtype=int) #o tamanho de todas as sequencias eh igual
											  #A
											  #C
											  #G
											  #T

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

	for baseVector in positionCountMatrix:
		for base in baseVector:
			print base,
		print '\n'


def finalMotif(vector):
	return
def similarity(vector):
	return
def support(vector):
	return
def complexity(vector):
	return
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
	
	pcm = positionCountMatrix(dna_sequences)
	return
# ---- SOLVE TEST CASE WITH ARTIFICIAL BEE COLONY ALGORITHM

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
	



	readFasta("assets/Real/dm01r.fasta")
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
