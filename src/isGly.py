def is_PNGS(seq, site):
	#n-linked glycosylation motif: NXS/T/C
	triple = seq[site:site+3]
	
	#check if this substitution can cause glycosylation
	if (triple[0] == 'N') and (triple[2] == 'S' or triple[2] == 'T' or triple[2] == 'C'):
		return 1 

	#not PNGS
	return 0