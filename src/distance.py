def is_PNGS(seq, site):
	#n-linked glycosylation motif: NXS/T/C
	triple = seq[site:site+3]
	
	#check if this substitution can cause glycosylation
	if (triple[0] == 'N') and (triple[2] == 'S' or triple[2] == 'T' or triple[2] == 'C'):
		return 1 

	#not PNGS
	return 0
	
def calc_similarities(byYear, firstY, lastY, site, vac):
	Sims = []
	
	for y in range(len(byYear)):
		Sim = 0
		for seq in byYear[y]:
			s = 0
			for l in site:
				if seq[l-1] == vac[l-1]:
					s += 1
			Sim += s
		try:
			Sim = 1.0*Sim/len(byYear[y])
			Sims.append(Sim)
		except ZeroDivisionError:
			Sims.append("NA")
	return Sims