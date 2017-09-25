def read_epitope_shih(epifName):
	sites = ['a','b','c','d']
	epitopes = [[] for i in range(len(sites)) ]
	
	epif = open(epifName, "r")
	for line in epif:
		each = line.split("\n")[0].split("\t")
		idx = sites.index(each[1])
		epitopes[idx].append(int(each[0]))
			
	epif.close()
	
	return epitopes










	
	
	
	
	
	
	