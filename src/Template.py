from ModObj import ModObj
from Alignment import Alignment
from PDB import PDB

class Template(ModObj):


	def __init__(self,id=None,fasta=None,sequence=None,length=0,alignments=[],pdb=None):
		self.id = id
		self.fasta = fasta
		self.sequence = sequence
		self.length = length
		self.alignments = alignments
		self.pdb = pdb