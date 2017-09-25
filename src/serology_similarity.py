import os
import pandas

###########################################################################
# Take HIVE1314 data, aggregate same birth year, get average pre and post Texas titer

firstY = 1937
lastY = 1999

###	read hive data
hivefName = os.path.normpath("../data/HIVE1314.CSV")
hiveoutfName = os.path.normpath("../data/HIVE1314_import.csv")

hive = pandas.read_csv(hivefName)

hive_out = open(hiveoutfName, "w")
hive_out.write("id,pre,post\n")

###	put titers by year of birth
byYOB_pre = [[] for i in range(firstY, lastY+1)]
byYOB_post = [[] for i in range(firstY, lastY+1)]

for i in range(len(hive['id'].values)):
	id = str(hive['id'].values[i])
	try:
		
		ex1h3pre = float(hive['ex1-h3-pre'].values[i])
		ex2h3pre = float(hive['ex2-h3-pre'].values[i])
		ex1h3post = float(hive['ex1-h3-post'].values[i])
		ex2h3post = float(hive['ex2-h3-post'].values[i])
		
		if type(ex2h3post) != type(1.0):
			print (ex2-h3-post)
		
		pre = 1.0*(float(hive['ex1-h3-pre'].values[i]) + float(hive['ex2-h3-pre'].values[i])) / 2
		post = 1.0*(float(hive['ex1-h3-post'].values[i]) + float(hive['ex2-h3-post'].values[i])) / 2
		yob = int(hive['yob'].values[i])
		byYOB_pre[yob-firstY].append(pre)
		byYOB_post[yob-firstY].append(post)
	except ValueError:
		continue
	
###	get average pre and post titers of each cohort
for y in range(len(byYOB_pre)):
	try:
		avg_pre = 1.0*sum(byYOB_pre[y])/len(byYOB_pre[y])
		avg_post = 1.0*sum(byYOB_post[y])/len(byYOB_post[y])
		oneline = str(y+firstY)+","+str(avg_pre)+","+str(avg_post)
	except ZeroDivisionError:
		oneline = str(y+firstY)+",NA,NA"
	hive_out.write(oneline+"\n")
	
hive_out.close()

################################################################################
# Calculate A similarity, B similarty, #gly at A, #gly at B, gly similarity, HA2 similarity
import epi_similarity
import distance

###	read vaccine strain, read H3N2 sequences and put them in list by year
vacfName = os.path.normpath("../data/Texas2012_aa.fas")
vacf = open(vacfName, "rU")
for line in vacf:
	if line.find(">") > 0:
		continue
	else:
		vac = line.split('\n')[0]
		
byYear = [[] for i in range(firstY, lastY+1)]
seqfName = os.path.normpath("../data/ncbi_aligned_6812_AA.fas")
seqf = open(seqfName, "rU")
for line in seqf:
	if line.find(">") >= 0:
		each = line.split("\n")[0].split("|")
		y = int(each[2].split("/")[0])
	else:
		if y < firstY or y > lastY:
			continue
		byYear[y-firstY].append(line.split("\n")[0])
			
############todo: remove ambiguous		
		
###	determine epitope A and B by Shih et al. 2009 PNAS

shihfName = os.path.normpath("../data/shih_epitope.txt")	
epitopes = epi_similarity.read_epitope_shih(shihfName)
epitopeA = epitopes[0]
epitopeB = epitopes[1]

###	determine HA2
HA2 = range(330,550+1)

### calc epitope similarities
A_similarities = distance.calc_similarities(byYear, firstY, lastY, epitopeA, vac)
B_similarities = distance.calc_similarities(byYear, firstY, lastY, epitopeB, vac)

### calc HA2 similarities
HA2_similarities = distance.calc_similarities(byYear, firstY, lastY, HA2, vac)

### calc glycosylation similarities


### write sero-similarities file

sero_simfName = os.path.normpath("../data/HIVE1314_similarities.csv")
sero_simf = open(sero_simfName, "w")
sero_simf.write("id,pre,post,Asimilarity,Bsimilarity,HA2similarity\n")

hive_out = open(hiveoutfName, "rU")

sero = []
sero_sim = ['' for i in range(firstY, lastY)]

for line in hive_out:
	if line.find("id") >= 0: 
		continue
	else:
		sero.append(line.split("\n")[0])

for y in range(0, lastY-firstY+1):
	oneline = sero[y]+","
	oneline += str(A_similarities[y])+"," +str(B_similarities[y])
	oneline += ","+str(HA2_similarities[y])
	sero_simf.write(oneline+"\n")
	




























