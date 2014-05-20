def convertToGFF3 (filename,outputName,seqname): 
	out = open(outputName,'w')

	with open(filename, "r") as f:
		contents = f.read()

	lines = contents.split("\n")
	#get first line which contains sequence information
	seqname = seqname + '\t'
	source = 'glimmer3' + '\t'
	feature = 'gene' + '\t'

	for i in xrange (0, len(lines)-1):
		#ignore all other headers
		if lines[i][0] == ">" or lines[i][0] == "#": 
			continue

		cols = lines[i].split(" ")
		cols = [i for i in cols if i != '']

		if len(cols) == 0:
			continue

		ID = 'ID=' + 'glimmer_' + cols[0]
		start = cols[1] + '\t'
		stop = cols[2] + '\t'
		strand = (cols[3])[0] + '\t'
		frame = str(int((cols[3])[1]) - 1) + '\t'
		if len(cols)>=5:
			score = cols[4]+ '\t'
		else:
			score = '.'+'\t'
		out.write(seqname + source + feature  + start + stop + score + strand + frame + ID +  "\n")
	
	out.close()

#convertToGFF3("ribotrimmedfiltered_rRNA_tRNA-PeaksTest.predict", "ribotrimmedfiltered_rRNA_tRNA-PeaksTest.gff3", "ecolik12")
