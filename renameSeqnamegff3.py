
class renameSeq(object):
	def __init__(self,filename,outputName,seqname):
		self.filename = filename
		self.seqname = seqname
		self.outputName = outputName

	def renameSeqnameGFF3 (self): 
		out = open(self.outputName,'w')
	
		with open(self.filename, "r") as f:
			contents = f.read()
	
		lines = contents.split("\n")
		#get first line which contains sequence information
		self.seqname = self.seqname + '\t'
	
		for i in range (1, len(lines)):
			#ignore all other headers
			if lines[0] == ">": 
				continue
	
			cols = lines[i].split("\t")
			cols = [i for i in cols if i != '']
			print cols
			if len(cols) < 9:
				continue
			data = '\t'.join(cols[1:])
	
			out.write(self.seqname + data + '\n')
		
		out.close()
	
