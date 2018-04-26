import csv
import subprocess

def analyze(patients, firstVar, secondVar, row, multiPats):
	grepCmnd1 = firstVar+'.*'+secondVar+'\s(.*?)\s(.*)\s(.*)$'
	grepCmnd2 = secondVar+'.*'+firstVar+'\s(.*?)\s(.*)\s(.*)$'
                                
	sep = '\t'

	patientSmartPhaseOut = './outfiles/'+patients+'.smartPhase.out'

	try:
		grep1 = subprocess.check_output(['grep', '-P', grepCmnd1, patientSmartPhaseOut])
		if(grep1 != ""):
			cols = grep1.split()
			outString = sep.join(row)
			outString = outString+"\t"+cols[2]+"\t"+cols[3]+"\t"+cols[4]+"\n"
			if(not multiPats):
				outFile.write(outString)
			else:
				return cols[2]+"|"+cols[3]+"|"+cols[4]
	except:
		try:
			grep2 = subprocess.check_output(['grep', '-P', grepCmnd2, patientSmartPhaseOut])
			if(grep2 != ""):
				cols = grep2.split()
				outString = sep.join(row)
				outString = outString+"\t"+cols[2]+"\t"+cols[3]+"\t"+cols[4]+"\n"
				if(not multiPats):
					outFile.write(outString)
				else:
					return cols[2]+"|"+cols[3]+"|"+cols[4]
		except:
			
			print "For patient "+patients+", couldn't find var: "+grepCmnd1
			pass


with open('NextSeqCohort_potCompHetSample.tsv', 'rb') as csvfile:
	linereader = csv.reader(csvfile, delimiter='\t', quotechar='"')
	header = linereader.next()
	outFile = open('NextSeqCohort_potCompHetSample.phased.tsv', 'w')
	sep = '\t'
	outStr = sep.join(header)
	outStr = outStr+"\tCompound Het. Prediction\tConfidence\n"
	outFile.write(outStr)
	for row in linereader:
		firstVar = row[1]
		secondVar = row[2]
		if(firstVar != "" and secondVar != ""):
			patients = row[0]
			if(',' in patients ):
				multiPatients = patients.split(",")
				predScore = ''
				predConf = ''
				predReads = ''
				for pat in multiPatients:
					pat = pat.strip()
					valsString = analyze(pat, firstVar, secondVar, row, True)
					if(valsString is None):
						continue
					valsSplit = valsString.split("|")
					predScore = predScore+valsSplit[0]+","
					predConf = predConf+valsSplit[1]+","
					predReads = predReads+valsSplit[2]+","
				predScore = predScore[:-1]
				predConf = predConf[:-1]
				predReads = predReads[:-1]
				#predScore = predScore+'"'
				#predConf = predConf+'"'
				outString = sep.join(row)
				outString = outString+"\t"+predScore+"\t"+predConf++"\t"+predReads+"\n"
				outFile.write(outString)
			else:
				analyze(patients, firstVar, secondVar, row, False)
		else:
			outString = sep.join(row)
			outString = outString+"\t\t\n"
			outFile.write(outString)
