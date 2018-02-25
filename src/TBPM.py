from Protein import Protein

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import re
import os
import shutil

targets = []

debug = True

skipBlanks = False
modellerCompatibilityMode = True

def main():
	for target in targets:
		
		print()
		print("Target:",target)
		print()
		target.getTemplates()
		target.getPDBs()
		target.alignPDB(skipBlanks = skipBlanks, modellerCompatibilityMode = modellerCompatibilityMode)

		# If we are skipping blanks we can use Scwrl4
		# Otherwise, fake the operation for now
		if skipBlanks:
			os.system("Scwrl4 -i targets/%s.pdb -o targets/%s_chain.pdb" % (target,target))
		else:
			shutil.copyfile('targets/%s.pdb' % (target), 'targets/%s_chain.pdb' % (target))

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Error: Expected 2+ arguments, received {0}".format(len(sys.argv)))
		sys.exit()
	i = 1
	while i < len(sys.argv):
		input = sys.argv[i]
		if ".fasta" in input:
			for seq_record in SeqIO.parse(input, "fasta"):
				seq_record.seq.alphabet = IUPAC.protein
				targets.append(Protein(seq_record.id,seq_record.seq,debug=debug))
		else:
			id = input
			i = i+1
			input = sys.argv[i]
			targets.append(Protein(id,Seq(input,IUPAC.protein),debug=debug))
		i = i+1

	main()