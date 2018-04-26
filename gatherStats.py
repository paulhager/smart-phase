import csv
import subprocess


with open('NextSeqCohort_potCompHetSample.phased.tsv', 'rb') as csvfile:
	linereader = csv.reader(csvfile, delimiter='\t', quotechar='"')
	header = linereader.next()
	sep = '\t'
	totalVars = 0
	totalCis = 0
	totalTrans = 0
	totalInnoc = 0
	weakVariant = 0
	toRemoveVar = 0
	totalDifference = 0
	withinOneRead = 0
	within500 = 0
	phased = 0
	phasedWithinOneRead = 0
	proximalUnphased = 0
	okCisCount = 0
	cisPhaseOneRead = 0
	for row in linereader:
		if("pa" in row[0]):
			firstVarSplit = row[1].split("-")
			secondVarSplit = row[2].split("-")
			distance = abs(int(secondVarSplit[1])-int(firstVarSplit[1]))
			totalDifference += distance
			if(distance < 128 and distance != 0 and row[3] == "2"):
				expectedReads = (150 - distance + 1)/3
				if(float(row[5]) > (expectedReads-1)):
					proximalUnphased += 1
					print sep.join(row)
			if(distance < 150):
				withinOneRead += 1
			if(distance < 500):
				within500 += 1
			totalVars += 1
			cis = row[3].count("0")
			if(cis == 1 and float(row[4]) >= 0.1):
				okCisCount += 1
			trans = row[3].count("1")
			innoc = row[3].count("4")
			innoc += row[3].count("5")
			innoc += row[3].count("6")
			totalInnoc += innoc
			totalCis += cis
			totalTrans += trans
			if(innoc + cis + trans > 0):
				phased += 1
				if(abs(int(secondVarSplit[1])-int(firstVarSplit[1])) < 150):
					phasedWithinOneRead += 1
					if cis == 1:
						cisPhaseOneRead += 1
	print "Total Vars: "+str(totalVars)
	print "Total Cis: "+str(totalCis)
	print "Total Trans: "+str(totalTrans)
	print "Total Innoc: "+str(totalInnoc)
	print "Cis > 0.1: "+str(okCisCount)
	print "Proximal Unphased: "+str(proximalUnphased)
	print "Average Distance: "+str(totalDifference/totalVars)
	print "Within average read: "+str(withinOneRead)
	print "Wtihin 500 bp: "+str(within500)
	print "Phased vars: "+str(phased)
	print "Phased within 1 read: "+str(phasedWithinOneRead)
	print "Cis Phase one read: "+str(cisPhaseOneRead)

