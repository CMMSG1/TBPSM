from ModObj import ModObj

class Atom(ModObj):

	def __init__(self,atomName='    ',altLoc=' ',resName='---',chainId=' ',
		codeForInsertion=' ',xcoord=0.0,ycoord=0.0,zcoord=0.0,
		occ=0.0,temp=0.0,elemSym='  ',charge='  ',missing=False):
		self.atomName = atomName
		self.altLoc = altLoc
		self.resName = resName
		self.chainId = chainId
		self.codeForInsertion = codeForInsertion
		self.xcoord = xcoord
		self.ycoord = ycoord
		self.zcoord = zcoord
		self.occ = occ
		self.temp = temp
		self.elemSym = elemSym
		self.charge = charge
		self.missing = missing