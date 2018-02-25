from ModObj import ModObj

class Alignment(ModObj):

	def __init__(self,id=None,title=None,expect=0,score=0,
		identities=0,similarity=0,target=None,targetstart=0,
		match=None,template=None,templatestart=0,length=0):
		self.id = id
		self.title = title
		self.expect = expect
		self.score = score
		self.identities = identities
		self.similarity = similarity
		self.target = target
		self.targetstart = targetstart
		self.match = match
		self.template = template
		self.templatestart = templatestart
		self.length = length