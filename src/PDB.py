
from ModObj import ModObj

class PDB(ModObj):


	def __init__(self,infos={},lines=[]):
		self.infos = infos
		self.lines = lines