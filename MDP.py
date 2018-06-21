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
from Hive import HiveMotif
from Hive import Utilities

random.seed(time.time()) 




def readFasta(filePath):
	
	dna_sequences = []
	fasta_sequences = SeqIO.parse(open(filePath),'fasta')
	i = 0
	for fasta in fasta_sequences:
		name,sequence = fasta.id,str(fasta.seq)
		dna_sequences.append(sequence)
		print(">Seq",i,dna_sequences[len(dna_sequences)-1])
		i += 1
	return dna_sequences


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
			print("start = ",motifStart)
			subSequence = sequence[motifStart:motifStart+motifSize]
			dna_subSequences.append(subSequence)
			print(subSequence)
			i += 1
	else:
		print("Sequencias de tamanho insuficiente:",sequenceSize,"<",lowerLimit)

	return dna_subSequences


def run():
    ###############readFasta(sys.argv[1]) #1 = path name
    dna_sequences = readFasta("assets/Real/mus01r.fasta")

    dnaLength = len(dna_sequences[0])
    model = HiveMotif.BeeHive(
    lower = [0],
    upper      = [64],
    dna_sequences = dna_sequences,
    fun        = None,
    numb_bees  = 50,
    max_itrs   = 100,
    max_trials = 10)
    cost = model.run()
    print("Fitness Value ABC: {0}".format(model.best))
    #Utilities.ConvergencePlot(cost)
"""
    solutionData = randomSubSequences(dna_sequences)
    solutionInstance = solutionData[0]
    dna_subSequences = solutionData[1]
    print("Solution Instance: ",solutionInstance)

    pcm = positionCountMatrix(dna_subSequences)
    consensus = consensusMotif(pcm)
    dna_approvedSequences = thresholdConsensus(dna_subSequences,consensus)
    finalPcm = positionCountMatrix(dna_approvedSequences)
    finalPfm = positionFrequencyMatrix(finalPcm)
    print("Motif:")
    finalMotif = consensusMotif(finalPcm)
    motifSimilarity = similarity(finalPfm)
    motifComplexity = complexity(finalMotif)
    motifSupport= len(dna_approvedSequences)
    biased = isBiased(motifSupport,len(dna_subSequences))
    print("Similarity:",motifSimilarity)
    print("Complexity:",motifComplexity)
    print("Support:",motifSupport)
    print("Biased?",biased)
"""

if __name__ == "__main__":
	run()


# ---- END
