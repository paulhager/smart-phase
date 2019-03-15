#  Copyright (C) 2018 the SmartPhase contributors.
#  Website: https://github.com/paulhager/smart-phase
#
#  This file is part of the SmartPhase phasing tool.
#
#  The SmartPhase phasing tool is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import vcf
import csv

goldenVCFPath = sys.argv[1]
phaseCallsTsv = sys.argv[2]
sampleName = sys.argv[3]

goldenVCFReader = vcf.Reader(filename = goldenVCFPath)
with open(phaseCallsTsv) as tsvFile:
  linereader = csv.reader(tsvFile, delimiter='\t')
  header = next(linereader)
  print ("\t".join(header)+"\tActual")
  for row in linereader:
    firstVar = row[1]
    secondVar = row[2]
    call = row[3]
    if(call != '4'):
      firstVarSplit = firstVar.split("-")
      secondVarSplit = secondVar.split("-")
      contig = firstVarSplit[0]
      firstVarStart = firstVarSplit[1]
      firstVarRef = firstVarSplit[2]
      firstVarAlt = firstVarSplit[3]
      secondVarStart = secondVarSplit[1]
      secondVarRef = secondVarSplit[2]
      secondVarAlt = secondVarSplit[3]
      goldenVariants = goldenVCFReader.fetch(contig, start=int(firstVarStart)-1, end=int(secondVarStart)+1)
      firstVarGoldenGT = ""
      secondVarGoldenGT = ""
      for variant in goldenVariants:
        if(str(variant.POS) == firstVarStart and variant.REF == firstVarRef and variant.ALT[0] == firstVarAlt):
          firstVarGoldenGT = variant.genotype(sampleName)
        elif(str(variant.POS) == secondVarStart and variant.REF == secondVarRef and variant.ALT[0] == secondVarAlt):
          secondVarGoldenGT = variant.genotype(sampleName)
      if(firstVarGoldenGT and secondVarGoldenGT):
        if(firstVarGoldenGT.gt_nums == secondVarGoldenGT.gt_nums):
          goldenCall = '1'
        else:
          goldenCall = '2'
        if(goldenCall != call):
          print ("\t".join(row)+"\t"+str(goldenCall))
      else:
        print ("Error! Couldn't find both vars in golden file! Vars were: "+firstVar+" || "+secondVar)
