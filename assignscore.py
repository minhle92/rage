def AssignScore (filename,scores): 
	
	with open(filename, "r") as f:
		contents = f.read()

	lines = contents.split("\n")
	#get first line which contains sequence information
	out = open(filename,'w')
	for i in xrange (0, len(lines)-1):
		#ignore all other headers
		if lines[i][0] == ">" or lines[i][0] == "#": 
			continue
		cols = lines[i].strip()
		cols = cols.split("\t")
		cols[5] = scores
		
		data = '\t'.join(cols)
		out.write(data + '\n')
	
	out.close()
	
def writeResults2gff3(results, outFile, genomeinfo,source):
	out = open(outFile,'w')
	for index in xrange(len(results)):
		if results[index][0]>results[index][1]:
			signs = "-"
		else:
			signs = "+"
		results[index][0]= str(results[index][0])
		results[index][1]= str(results[index][1])
		results[index][2]= str(results[index][2])
		data = [genomeinfo,source,"gene"]+results[index]+[signs,".","ID=Combiner_orf"+str(index)]
		data = '\t'.join(data)
		out.write(data + '\n')
	out.close()

#AssignScore("rnapred.gff3" , '3')
#writeResults2gff3([[3,4,2.2],[2,6,3.3]],"testing.gff3","abc","test")
