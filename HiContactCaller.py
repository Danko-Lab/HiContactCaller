import sys
import numpy as np
import statsmodels.stats.multitest as mt
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector as ivect
import rpy2.robjects as robjects
import rpy2
import gzip

mass = importr('MASS')
stats = importr('stats')
base = importr('base')
    
peakFile = sys.argv[1] # File of "bait" peaks to analyze (format: Chromosome <\t> bait-center-position)
preyFile =  sys.argv[2] # File of "prey" PREYs to analyze (format: Chromosome <\t> prey-center-position)
conFile = sys.argv[3] # Contact file in juicer short format (given from merged_nodup.txt)
out = sys.argv[4] # Path to output file
DIST = int(sys.argv[5]) # Half size of "prey" search window around "bait" position
CAP = int(sys.argv[6]) # Half size of "bait" and "prey" contact capture window
DRP = int(sys.argv[7]) # Length of drop in contact frequency after CAP - probably a function of fragment size ditribution and allowed intra-fragment contacts

loci, peaks, PREYs, contactProbabilities, pvalues = [],[],[],[],[]

print "Getting loci for contact analysis"

#The following loop fills up the list of loci surrounding "bait" peaks and adds empty list to their corresponding 
#positions in the list of contacts and "preys" to list those contacts and PREYs that falls within these loci
for locus in open(peakFile): # locus line format: Chromosome <\t> peak-center-position>
    locus = locus.strip().split()
    start = max(0, int(locus[1]) - DIST)
    stop = int(locus[1]) + DIST
    loci.append((locus[0],start,stop,int(locus[1])))
    PREYs.append([])
    

print "Sorting Preys"

#The following loop assign "preys" to loci surrounding "bait" peaks so that each "prey"
#is writen as its distance from the center of the "bait" peak, the chromosome and its position
ct = 0
for prey in open(preyFile): # prey line format: Chromosome <\t> prey-center-position>
    split = prey.strip().split()
    chrom, pos = split[0], int(split[1])
    for i in range(len(loci)):
        if chrom == loci[i][0]:
            if loci[i][1] < pos < loci[i][2]:
                distance = pos - loci[i][3]
                PREYs[i].append((distance,chrom,pos))
    ct += 1
    if i%100 == 0: 
        print ct

print out[-6:-4], "Scanning",len(loci),"loci for interactions!"

# Calculate observed and expected distributions for each locus
print 'Analyzing', len(loci), 'loci'
counter = 0
for locus in open(peakFile): # locus line format: <chr>\t<position = center of peak>
    run = False
    expect, expect_close, plus, minus = [],[],[],[]
    locus = locus.strip().split()
    chrom, position = locus[0], int(locus[1])
    chrom = chrom.strip("chr")
    
    # read the contacts around loci from the splited temp_contact_file
    try:
        with gzip.open(conFile+"_IN_"+peakFile+"_folder/temp_"+str(counter)+".gz") as contact_fp:
        #with open(conFile+"_IN_"+peakFile+"_folder/temp_"+str(counter)) as contact_fp:
            contacts= np.array([l.strip().split() for l in contact_fp.readlines()])
            print "contacts with locus number", contacts[0][0]
            contacts = contacts[:,1:4]
            # Calculate parameters for expected distribution
            distance = contacts[:,2].astype(int) - contacts[:,1].astype(int)
            expect += list(distance)
            expect_close += [x for x in distance if x <= CAP + DRP]
            run = True
    except IOError:
        expect.append(1)
        expect_close.append(1)
    empDist = ivect(expect)
    robjects.r.assign("empDist",empDist)

    if len(expect_close) > 0:
        empDist_close = ivect(expect_close)
        robjects.r.assign("empDist_close",empDist_close)
        PC = robjects.r("ecdf(empDist_close)") # Empirical Cumulative Distribution Function for all contacts' distances in the locus within the capture window + drop length
        robjects.r.assign("PC",PC)

    P = robjects.r("ecdf(empDist)") # Empirical Cumulative Distribution Function for all contacts' distances in the locus
    robjects.r.assign("P",P)

    baitStart, baitStop = (position-CAP), (position+CAP)
    # Get all plus and minus diractions distances of interactions with positions within a 5kb (CAP*2) window from (+ and -) the peak's center
    if run:
		for contact in contacts:
			if baitStart <= int(contact[1]) <= baitStop:
				distance = int(contact[2])-int(contact[1])
				plus.append(distance)
			if baitStart <= int(contact[2]) <= baitStop:
				distance = int(contact[2])-int(contact[1])
				minus.append(distance)

		# test plus strand contacts with promoters
		for prey in PREYs[counter]:
			if prey[0] > 0: #making sure to exclude self-interactions
				preyChrom, preyPosition, preyDist = prey[1],prey[2], prey[0]
				preyStart, preyStop = (preyDist - (CAP)), (preyDist + (CAP))

				if (preyDist - (CAP + DRP)) > 0:
					expStart, expStop = max(CAP + DRP + 1,(preyDist - (CAP))), (preyDist + (CAP)) + max(0, ((CAP + DRP + 1) - (preyDist - CAP))))
					robjects.r.assign("expStart",expStart)
					robjects.r.assign("expStop", expStop)
					z = robjects.r("seq(expStart,expStop,by=1)")
					robjects.r.assign("z",z)
					p = robjects.r("P(z)")
					exp = p[-1]-p[0]
					expected = exp * len(plus)

				else:
					if preyDist > CAP:
						expStartC, expStopC = preyDist - CAP, CAP + DRP
					else:
						expStartC, expStopC = 0, CAP
					robjects.r.assign("expStartC",expStartC)
					robjects.r.assign("expStopC", expStopC)
					zC = robjects.r("seq(expStartC,expStopC,by=1)")
					robjects.r.assign("zC",zC)
					pc = robjects.r("PC(zC)")
					expC = pc[-1]-pc[0]
					if preyDist > CAP:
						expStart, expStop = CAP + DRP + 1, CAP + DRP + preyDist
					else:
						expStart, expStop = CAP + 1, CAP + preyDist
					robjects.r.assign("expStart",expStart)
					robjects.r.assign("expStop", expStop)
					z = robjects.r("seq(expStart,expStop,by=1)")
					robjects.r.assign("z",z)
					p = robjects.r("P(z)")
					exp = p[-1]-p[0]
					expected = (expC * len([x for x in plus if x <= CAP + DRP])) + (exp * len(plus))

				observed = 0
				for distance in plus:
					if preyStart <= distance <= preyStop:
						observed += 1

				v = robjects.FloatVector([observed, expected, (len(plus)-observed), (len(plus)-expected)])
				robjects.r.assign("v", v)
				test = robjects.r("matrix(v, nrow = 2)")
				robjects.r.assign("matrix", test)
				try:
					robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
					result = robjects.r("exact")
					p_val = result[0][0]
				except:
					p_val = 999.99
				contactProbabilities.append([chrom, position, preyChrom, preyPosition, preyDist, p_val, observed, expected])
				pvalues.append(p_val)

		# test minus strand contacts
		for prey in PREYs[counter]:
			if prey[0] < 0:
				preyChrom, preyPosition, preyDist = prey[1],prey[2], prey[0]
				preyStart, preyStop = (preyDist - (CAP)), (preyDist + (CAP))

				if (abs(preyDist) - (CAP + DRP)) > 0:
					expStart, expStop = max(CAP + DRP + 1, (abs(preyDist) - (CAP))), (abs(preyDist) + (CAP)) + max(0,((CAP + DRP + 1) - (abs(preyDist) - CAP))))
					robjects.r.assign("expStart",expStart)
					robjects.r.assign("expStop", expStop)
					z = robjects.r("seq(expStart,expStop,by=1)")
					robjects.r.assign("z",z)
					p = robjects.r("P(z)")
					exp = p[-1]-p[0]
					expected = exp * len(minus)

				else:
					if abs(preyDist) > CAP:
						expStartC, expStopC = abs(preyDist) - CAP, CAP + DRP
					else:
						expStartC, expStopC = 0, CAP
					robjects.r.assign("expStartC",expStartC)
					robjects.r.assign("expStopC", expStopC)
					zC = robjects.r("seq(expStartC,expStopC,by=1)")
					robjects.r.assign("zC",zC)
					pc = robjects.r("PC(zC)")
					expC = pc[-1]-pc[0]
					if abs(preyDist) > CAP:
						expStart, expStop = CAP + DRP + 1, CAP + DRP + abs(preyDist)
					else:
						expStart, expStop = CAP + 1, CAP + abs(preyDist)
					robjects.r.assign("expStart",expStart)
					robjects.r.assign("expStop", expStop)
					z = robjects.r("seq(expStart,expStop,by=1)")
					robjects.r.assign("z",z)
					p = robjects.r("P(z)")
					exp = p[-1]-p[0]
					expected = (expC * len([x for x in minus if x <= CAP + DRP])) + (exp * len(minus))


				observed = 0
				for distance in minus:
					if preyStart <= (-1*distance) <= preyStop:
						observed += 1

				v = robjects.FloatVector([observed, expected, (len(minus)-observed), (len(minus)-expected)])
				robjects.r.assign("v", v)
				test = robjects.r("matrix(v, nrow = 2)")
				robjects.r.assign("matrix", test)
				try:
					robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
					result = robjects.r("exact")
					p_val = result[0][0]
				except:
					p_val = 999.99

				contactProbabilities.append([chrom, position, preyChrom, preyPosition, preyDist, p_val, observed, expected])
				pvalues.append(p_val)

	counter += 1

out = open(sys.argv[4], "w")

corrected = mt.multipletests(pvalues, method='fdr_bh')
print 'writing'
p_count = 0
for prob in contactProbabilities:
    # output file format - chr <\t> bait position <\t> prey position <\t> p_val <\t> observed <\t> expected <\t> FDR
    outStr = str('chr' + str(prob[0])) + "\t" + str(prob[1]) +"\t"+ str(prob[3]) +"\t"+ str(prob[4]) +"\t"+ str(prob[5]) +"\t"+ str(prob[6]) +"\t"+ str(prob[7]) +"\t"+ str(corrected[1][p_count]) +"\n"
    out.write(outStr)
    p_count += 1

