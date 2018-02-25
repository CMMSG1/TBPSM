from Template import Template
from Alignment import Alignment
from PDB import PDB
from Atom import Atom

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "abcd@example.com"
from Bio import File

import sqlite3
import urllib
import io
import os
import time
import re
import sys
import json
import copy

rcsb_url = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s' 

amino3to1 = {
	"ALA":"A",
	"ARG":"R",
	"ASN":"N",
	"ASP":"D",
	"CYS":"C",
	"GLU":"E",
	"GLN":"Q",
	"GLY":"G",
	"HIS":"H",
	"ILE":"I",
	"LEU":"L",
	"LYS":"K",
	"MET":"M",
	"PHE":"F",
	"PRO":"P",
	"SER":"S",
	"THR":"T",
	"TRP":"W",
	"TYR":"Y",
	"VAL":"V",
	"MYL":"X",
	"---":"-"
}
amino1to3 = {
	"A":"ALA",
	"R":"ARG",
	"N":"ASN",
	"D":"ASP",
	"C":"CYS",
	"E":"GLU",
	"Q":"GLN",
	"G":"GLY",
	"H":"HIS",
	"I":"ILE",
	"L":"LEU",
	"K":"LYS",
	"M":"MET",
	"F":"PHE",
	"P":"PRO",
	"S":"SER",
	"T":"THR",
	"W":"TRP",
	"Y":"TYR",
	"V":"VAL",
	"X":"MYL",
	"-":"---"
}

class Protein:


	NMR_RESOLUTION = 3.5

	debug = False

	def __init__(self,pid=None,sequence=None,debug=False):
		self.pid = pid if pid else str(datetime.now())
		self.seq = sequence
		self.debug = debug
		self.templates = {}
		self.fastas = {}
		self.alignments = []
		self.pdb = None
		self.pdbs = {}
		
		self.fastasfolder = 'fastas-%s' % str(self.pid)
		self.alignmentsfolder = 'alignments-%s' % str(self.pid)
		self.templatesfolder = 'templates-%s' % str(self.pid)
		self.targetsfolder = 'targets'
		if not os.path.exists(self.fastasfolder):
			os.mkdir(self.fastasfolder)
		if not os.path.exists(self.alignmentsfolder):
			os.mkdir(self.alignmentsfolder)
		if not os.path.exists(self.templatesfolder):
			os.mkdir(self.templatesfolder)
		if not os.path.exists(self.targetsfolder):
			os.mkdir(self.targetsfolder)
		
		self.ex_resolution = re.compile(b'REMARK\s*2 RESOLUTION\.\s*([0-9\.]+|NOT APPLICABLE)\s+.*')

	def saveToDb(self):
		print("Save to db")

	def getTemplates(self):

		if self.debug:
			print(self.seq)

		result_handle = NCBIWWW.qblast("blastp","pdb",str(self.seq),expect=0.01)		
		blast_records = NCBIXML.parse(result_handle)
		if self.debug:
			print("BLAST Request Finished")
			print()
		
		for record in blast_records:			
			for alignment in record.alignments:				
				id = alignment.accession
				fasta = self.getFastaFromId(id)
				title = alignment.title
				length = alignment.length
				
				template = Template(
					id=id,fasta=fasta,sequence=title,
					length=length,alignments=[]
				)
				
				self.templates[id] = template				
				self.fastas[id] = fasta				
				for hsp in alignment.hsps:
					
					a = Alignment(
						id=id,title=title,expect=hsp.expect,score=hsp.score,
						identities=hsp.identities,similarity=(100*hsp.identities/len(self.seq)),
						target=hsp.query,targetstart=hsp.query_start,match=hsp.match,
						template=hsp.sbjct,templatestart=hsp.sbjct_start,length=length
					)
					
					targetfront = str(self.seq[:a.targetstart-1])
					targetend = str(self.seq[(a.targetstart+a.length):])
					a.target = ''.join(targetfront) + a.target + ''.join(targetend)
					a.length = len(a.target)
					
					templatefront = ['-']*(a.targetstart-1)
					templateend = ['-']*(len(self.seq)-(a.targetstart+a.length))
					a.template = ''.join(templatefront) + a.template + ''.join(templateend)

					self.templates[id].alignments.append(a)
					self.alignments.append(a)

		for id,fasta in self.fastas.items():
			fname = '%s/%s.fasta' % (self.fastasfolder,id)
			if not os.path.exists(fname):
				f = open(fname,'w')
				SeqIO.write(fasta,f,'fasta')
				f.close()

		for i,a in enumerate(self.alignments):
			fname = '%s/%s-%s.alignment' % (self.alignmentsfolder,a.id,str(i))
			if not os.path.exists(fname):
				f = open(fname,'w')
				json.dump(a.toJSON(),f)
				f.close()

		return self.templates.keys()

	def getFastaFromId(self,id):
		handle = Entrez.efetch(db="protein", rettype="fasta", id=id)
		frecord = SeqIO.read(handle, "fasta")
		frecord.id = str(id)		
		handle.close()
		return frecord

	def getPDBs(self):
		results = []
		for id in self.templates.keys():
			handle = self.getRemotePDBHandle(id.split('_')[0])
			lines,infos = self.parsePdbFromHandle(handle)
			pdb = PDB(infos=infos,lines=lines)
			if id in self.templates:
				self.templates[id].fasta.__dict__.update(infos)
			self.templates[id].pdb = pdb
			self.pdbs[id] = pdb
			fname = '%s/%s.pdb' % (self.templatesfolder,id)
			if not os.path.exists(fname):
				f = open(fname,'wb',1)
				f.writelines(lines)
				f.close()
			results.append(fname)
		return results

	def getRemotePDBHandle( self, id ):

		handle = urllib.request.urlopen( rcsb_url% (id) )
		uhandle = File.UndoHandle(handle)
		if not uhandle.peekline():
			raise BaseException( "Couldn't retrieve ", rcsb_url )
		return uhandle

	def parsePdbFromHandle(self, handle, first_model_only=True ):

		lines = []
		res_match = None
		infos = {}

		if type( handle ) is str:
			if len(handle) < 5000:
				raise BaseException( "Couldn't extract PDB Info." )
			handle =  io.StringIO( handle )

		for l in handle:
			lines += [ l ]

			res_match = res_match or self.ex_resolution.search( l )

			if first_model_only and l[:6] == b'ENDMDL':
				break

		if res_match:
			if res_match.groups()[0] == 'NOT APPLICABLE':
				infos['resolution'] = self.NMR_RESOLUTION
			else:
				infos['resolution'] = float( res_match.groups()[0] )

		return lines, infos

	def alignPDB(self, skipBlanks = True, modellerCompatibilityMode=False):
		
		alignments = self.alignments
		alignments.sort(key = lambda a: a.score,reverse=True)
		target = None
		template = None
		alignment = None
	
		for a in alignments:
			if a.similarity > 90: 
				continue
			else:
				alignment = a
				target = a.target
				template = a.template
				pdb = self.pdbs[a.id]
				break

		if self.debug:
			print("Selected alignment:",alignment.id)
			print()
			print("Selected template:",template)
			print()

		seq = []
		curres = '---'
		res = []
		first = True
		atomnum = 5
		resnum = 0
		for line in pdb.lines:
			
			line = line.decode()
			
			if line[:6] == 'ATOM  ':
				
				atomname = line[12:16]

				if atomname not in (' N  ',' CA ',' C  ',' O  ') or line[16] not in (' ','A'):
					continue

				if line[17:20] != curres or atomnum > 4:

					if atomnum < 5:
						if atomnum == 1:
							atom = Atom(atomName=' N  ',resName=curres,elemSym=' N',missing=True)
							res.append(atom)
							atomnum += 1
						if atomnum == 2:
							atom = Atom(atomName=' CA ',resName=curres,elemSym=' C',missing=True)
							res.append(atom)
							atomnum += 1
						if atomnum == 3:
							atom = Atom(atomName=' C  ',resName=curres,elemSym=' C',missing=True)
							res.append(atom)
							atomnum += 1
						if atomnum == 4:
							atom = Atom(atomName=' O  ',resName=curres,elemSym=' O',missing=True)
							res.append(atom)
							atomnum += 1
					if res: seq.append((amino3to1[curres],res))

					if first:
						first = False

						for k in range(1,int(line[22:26])):
							resname = '---'
							res = []
							res.append(Atom(atomName=' N  ',resName=resname,elemSym=' N',missing=True))
							res.append(Atom(atomName=' CA ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' C  ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' O  ',resName=resname,elemSym=' O',missing=True))
							seq.append((amino3to1[resname],res))
						resnum = int(line[22:26])-1

					resinc = int(line[22:26])-resnum
					if resinc != 1:

						for k in range(1,resinc):
							resname = '---'
							res = []
							res.append(Atom(atomName=' N  ',resName=resname,elemSym=' N',missing=True))
							res.append(Atom(atomName=' CA ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' C  ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' O  ',resName=resname,elemSym=' O',missing=True))
							seq.append((amino3to1[resname],res))

					res = []
					resnum = int(line[22:26])
					atomnum = 1

				curres = line[17:20]

				if atomnum == 1:
					if atomname != ' N  ':
						atom = Atom(atomName=' N  ',resName=curres,elemSym=' N',missing=True)
						res.append(atom)
						atomnum += 1
				if atomnum == 2:
					if atomname != ' CA ':
						atom = Atom(atomName=' CA ',resName=curres,elemSym=' C',missing=True)
						res.append(atom)
						atomnum += 1
				if atomnum == 3:
					if atomname != ' C  ':
						atom = Atom(atomName=' C  ',resName=curres,elemSym=' C',missing=True)
						res.append(atom)
						atomnum += 1
				if atomnum == 4:
					if atomname != ' O  ':
						atom = Atom(atomName=' O  ',resName=curres,elemSym=' O',missing=True)
						res.append(atom)
						atomnum += 1
						continue
				if atomnum > 4:
					continue
				atom = Atom(
					atomName=atomname,
					altLoc=line[16],
					resName=line[17:20],
					chainId=line[21],
					codeForInsertion=line[26],
					xcoord=float(line[30:38]),
					ycoord=float(line[38:46]),
					zcoord=float(line[46:54]),
					occ=float(line[54:60]),
					temp=float(line[60:66]),
					elemSym=line[76:78],
					charge=line[78:80]
				)
				res.append(atom)
				atomnum += 1

		seq.append((amino3to1[curres],res))


		if self.debug:
			output = []
			for s in seq:
				output.append(s[0])

		i = 0 
		j = 0 
		targetseq = [] # final sequence
		if (len(template)<len(target)):
			length=len(template)
		else:
			length=len(target)

		while i < length and j < len(seq):

			if '-' == template[i]:
				

				resname = amino1to3[target[i]]
				res = []
				res.append(Atom(atomName=' N  ',resName=resname,elemSym=' N',missing=True))
				res.append(Atom(atomName=' CA ',resName=resname,elemSym=' C',missing=True))
				res.append(Atom(atomName=' C  ',resName=resname,elemSym=' C',missing=True))
				res.append(Atom(atomName=' O  ',resName=resname,elemSym=' O',missing=True))
				targetseq.append(res)
				i += 1
				if '-' == seq[j][0]:
					j += 1
				continue
			
			if '-' == target[i]:
				i += 1
				j += 1
				continue
			
			resname = amino1to3[target[i]]
			res = []
			
			for atom in seq[j][1]:
				if not atom.missing:
					a = copy.copy(atom)
					a.resName = resname
					res.append(a)
			
			targetseq.append(res)
			i += 1
			j += 1

		flag=False
		for atom in targetseq:
			if atom != []:
				 flag=True 
				 break

		if flag:
			lines = []
			lines.append("REMARK Template for target %s\n" % (self.pid))
			residues_needed_for_loop_modelling = []
			residues_needed_for_loop_modelling.append("Loops needed for target %s\n" % (self.pid))
        
			i_atom = 0
			i_residue = 0
			last_residue_for_loop_modelling = 0
			for res in targetseq:
				i_residue += 1
				for atom in res:
					i_atom += 1
					if atom.missing:
						if (last_residue_for_loop_modelling!=i_residue):
							residues_needed_for_loop_modelling.append("%i, " % i_residue)
							last_residue_for_loop_modelling = i_residue
						
						if (modellerCompatibilityMode):
							atom.xcoord=(i_atom%4)*2.0
					
					if not skipBlanks or not atom.missing:
						lines.append('ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (
							i_atom,atom.atomName,atom.altLoc,atom.resName,' ',i_residue,atom.codeForInsertion,
							atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.temp,atom.elemSym,atom.charge
						))
					
			fname = '%s/%s.pdb' % (self.targetsfolder,self.pid)
			f = open(fname,'w')
			f.writelines(lines)
			f.close()

			print("Remark Finish!")

		else:
			print("Match error!We can't get a PDB file")


	def __str__(self):
		return str(self.pid)
