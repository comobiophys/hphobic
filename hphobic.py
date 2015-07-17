from numpy import *
from prody import *
from os import listdir
import glob
import os
import sys
import os.path
from subprocess import check_output
from subprocess import call
import subprocess
import math
import getopt

AminoAcid = ["ala", "arg", "leu", "lys", "met", "gln" , "ile", "trp", "phe", "tyr", "cys" , "val", "asn" , "ser" , "his" , "glu", "thr" , "asp" , "gly" , "pro"]
Aminoacid = [ "a",   "r",   "l",   "k",   "m",   "q" ,   "i",   "w",   "f",   "y",   "c" ,   "v",   "n" ,   "s" ,   "h" ,   "e",   "t" ,   "d" ,   "g" ,   "p"]
HphoIndex = [  41 ,  -14 ,  97  , -23  ,  74  ,  -10  ,  99  ,  97  ,  100 ,   63 ,   49  ,  76  ,  -28  ,   -5  , 	  8  ,  -31 ,   13  ,  -55  ,    0  ,  -46 ] 

# Global variable
chain = ''
pdb = ''
rcsb = ''
output = ''

HPchain = []
HPratio = []
HPindex = []

def main(argv):
	global chain, pdb, rcsb, output

	try:
		opts, args = getopt.getopt(argv,"hc:p:r:o:",["chain=","pdb file=","rcsb=","output="])
	except getopt.GetoptError:
		print 'hphobic.py -c <Chain Sequence> -p <PDB File> -r <RCSB Id> -o <Output>'
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
		  print 'hphobic.py -c <Chain Sequence> -p <PDB File> -r <RCSB Id> -o <Output>'
		  sys.exit()
		elif opt in ("-c", "--chain"):
		  chain = arg
		elif opt in ("-p", "--pdb"):
		  pdb = arg
		elif opt in ("-r", "--rcsb"):
		  rcsb = arg
		elif opt in ("-o", "--output"):
		  output = arg

def isHyrophobic(aa):
	aa = aa.lower()

	if(aa == 'phe' or aa == 'f'):
		return True
	elif(aa == 'ile' or aa == 'i'):
		return True
	elif(aa == 'tpr' or aa == 'w'):
		return True
	elif(aa == 'leu' or aa == 'l'):
		return True
	elif(aa == 'val' or aa == 'v'):
		return True
	elif(aa == 'met' or aa == 'm'):
		return True
	elif(aa == 'tyr' or aa == 'y'):
		return True
	elif(aa == 'cys' or aa == 'c'):
		return True
	elif(aa == 'ala' or aa == 'a'):
		return True
	else: return False

def chainMethod():
	global HPratio, HPindex
	HPchain.append('main')
	HPratio.append(0.0)
	HPindex.append(0.0)

	for c in chain:
		if isHyrophobic(c) == True:
			HPratio[0] += 1
		HPindex[0] += HphoIndex[Aminoacid.index(c.lower())]

	# HPratio
	HPratio[0] /= len(chain)
	HPratio[0] *= 100
	print 'Hydrophobic ratio: '+str(HPratio[0])+'%'

	# HPindex
	HPindex[0] /= len(chain)
	print 'Hydrophobic index: '+str(HPindex[0])

def pdbMethod():
	# parse from pdb
	protein = parsePDB(pdb)
	calcByPrody(protein)

def rcsbMethod():
	# parse from rcsb server
	protein = parsePDB(rcsb)
	calcByPrody(protein)

def calcByPrody(protein):
	global HPratio, HPindex, HPchain
	listProtein = []

	hierView = protein.getHierView()
	chains = list(hierView)
	for y in range(0, len(chains)):
		chain = chains[y].select('protein')
		if chain is None: continue 
		#get hierview again
		chain = chain.getHierView()
		listProtein.append(list(chain)[0])

	# choose chains to calc
	print '\n*********************************'
	for x in range(0, len(listProtein)):
		print "{0:2d} => {1:10s} - {2:14s}{3:4s}{4:9s}".format(x, str(listProtein[x]), "protein with ", str(listProtein[x].numResidues()), ' residues')
	print 'Please select chains \nexample: \'chains: 0,1,2\' \n'
	char_list = raw_input("Enter chains: ")
	char_array = char_list.split(',')
	print '*********************************\n'

	# Calc hyrophobic
	for x in range(0, len(char_array)):
		if char_array[x].replace(" ", "") == "": continue
		index_x = int(char_array[x].replace(" ", ""))
		chain = listProtein[index_x]

		HPchain.append(str(listProtein[x]))
		HPratio.append(0.0)
		HPindex.append(0.0)

		for y in range(chain.getResnums()[0], chain.getResnums()[0]+chain.numResidues()):
			residue = chain.getResidue(y)
			resname = residue.getResname()
			if isHyrophobic(resname) == True:
				HPratio[x] += 1
			HPindex[x] += HphoIndex[AminoAcid.index(resname.lower())]

		print 'Chain '+HPchain[x]
		# HPratio
		HPratio[x] /= chain.numResidues()
		HPratio[x] *= 100
		print 'Hydrophobic ratio: '+ "{0:.2f}".format(HPratio[x]) +' %'

		# HPindex
		HPindex[x] /= chain.numResidues()
		print 'Hydrophobic index: '+ "{0:.2f}".format(HPindex[x])
		print '\n'


def outputStep():
	f = open(output,'w')
	for x in range(0, len(HPchain)):
		f.write('chain '+HPchain[x]+'\n')
		f.write('Hydrophobic ratio: '+"{0:.2f}".format(HPratio[x])+' %\n')
		f.write('Hydrophobic index: '+"{0:.2f}".format(HPindex[x])+'\n')
		f.write('\n')
	f.close()

if __name__ == "__main__":
	main(sys.argv[1:])

	if(chain != ''):
		chainMethod()

	if(pdb != ''):
		pdbMethod()

	if(rcsb != ''):
		rcsbMethod()

	if(output != ''):
		outputStep()
