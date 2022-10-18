#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### MODULES IMPORT ###
import argparse
import gzip
import collections
import datetime
import warnings
import pandas as pd
import numpy as np
from scipy.stats import linregress

### VARIABLES DEFINITION ###
GenTableFN = OutFN = MutaList = Host = DateBefore = DateAfter = Continent = Country = Region = Clade = ""
MutCutOff = DaysWindow = SeqLen = TopRank = 0
filters = []
    
### INITIALIZING TASKS ###
# Help about how to use the script
parser = argparse.ArgumentParser(description="02a_getGenotypeTable.py: Create a table containing genotype related data using genome's table as input (.tsv or .tsv.gz). \n Usage: ./02a_getGenotypeTable.py -i <GenomeTable(.tsv or .tsv.gz)> <option>\n", prefix_chars='-+')
parser.add_argument("-i", help="File name of the table containing all mutatated genomes and related information (.tsv or tsv.gz; required)")
parser.add_argument("-o", help="Output file name (tab separated)")
parser.add_argument("-dd", help="Number of days to consider in each time frame (default = 15)")
parser.add_argument("-mutN", help="Filer based on maximum number of mutations per genotype  (default = 100)")
parser.add_argument("-mutL", help="Filter mutations using a list of mutations separated by ',' (default = none; nucleotide, protein and annotations allowed)")
parser.add_argument("-host", help="Filter by host name (default = Human)")
parser.add_argument("-da", help="Filter by collection date after (YYYY-MM-DD) (default = none)")
parser.add_argument("-db", help="Filter by collection date before (YYYY-MM-DD) (default = none)")
parser.add_argument("-cont", help="Filter by continent name (default = none)")
parser.add_argument("-country", help="Filter by country name (default = none)")
parser.add_argument("-region", help="Filter by region name (default = none)")
parser.add_argument("-clade", help="Filter by Pangolin (Phylogenetic Assignment of Named Global Outbreak LINeages) clade (default = none)")
parser.add_argument("-seqlen", help="Filter by minimum sequence length (default = 29000)")
parser.add_argument("-nbases", help="Filter by the maximum number of degenerated bases found in sequence (default = 704)")
parser.add_argument("-toprank", help="Filter 'n' top ranked emerging genotypes (default = 0) (0 = no filter)")
parser.add_argument("-nocoldate", help="Do not output collection dates (True/False) (default = False)")
parser.add_argument("-nosubdate", help="Do not output submission dates (True/False) (default = False)")
parser.add_argument("-nogisaid", help="Do not output GISAID IDs (True/False) (default = False))")
parser.add_argument("-noclade", help="Do not output clade (True/False) (default = False)")
parser.add_argument("-noseqlen", help="Do not output sequence length (True/False) (default = False)")
parser.add_argument("-maxloc", help="Maximum number of locations to report (default = 50)")

#actually parse command line options 
args = parser.parse_args()
if args.i is not None:
	GenTableFN = args.i
else:
	print("  Please provide genomes in fasta format provided by GISAID (-i)")
	quit()

if args.o is not None:
	OutFN = args.o
	extension = args.o.rsplit(".")[-1]
	if extension != "gz" and  extension != "gzip":
	      OutFN = OutFN + '.gz'  
else:
        OutFN = "GenotypesTable.tsv.gz"

if args.mutN is not None:
	MutCutOff = int(args.mutN)
else:
	MutCutOff = 100

if args.mutL is not None:
	MutList = args.mutL.rsplit(",")
	print(MutList)
else:
	MutList = ""

if args.host is not None:
	RefGenome = args.host
else:
	Host = "Human"

if args.db is not None:
	DateBefore = args.db
else:
	DateBefore = ""
#	print("  Date before not defined so it is considered until the last date available. If you want to define an end date limite please use -db option)")

if args.da is not None:
	DateAfter = args.da
else:
	DateAfter = "2020-01-01" # imposes a starting date
#	print("  Date after not defined so it is considered from the first date available. If you want to define a starting date limite please use -da option)")

if args.dd is not None:
	DaysWindowDef = int(args.dd)
else:
	DaysWindowDef = 15
#	print("  The default value (15 days) will be applied for each time frame. If you want to change it use -dd option.")

DaysWindow = DaysWindowDef - 1 # correct the inclusion of the first and last day when creating time frames

if args.cont is not None:
	Continent = args.cont
else:
	Continent = ""
#	print("  Continent not defined. If you want to filter by continent, please use -cont option")

if args.country is not None:
	Country = args.country
	
else:
	Country = ""
#	print("  Country not defined. If you want to filter by country, please use -country option")

if args.region is not None:
	Region = args.region
else:
	Region = ""
#	print("  Region not defined. If you want to filter by region, please use -region option")
	
if args.clade is not None:
	Clade = args.clade
else:
	Clade = ""
#	print("  Clade not defined. If you want to filter by clade, please use -clade option")
	
if args.seqlen is not None:
	SeqLen = int(args.seqlen)
else:
	SeqLen = 29000
#	print("  Sequence length not defined. If you want to filter by sequence length, please use -seqlen option")

if args.nbases is not None:
	Nbases = int(args.nbases)
else:
	Nbases = 704
#	print("  Number of acceptable degenerated bases not defined. If you want to filter by the number of degenerated bases, please use -nbases option")


if args.toprank is not None:
	TopRank = int(args.toprank)
else:
	TopRank = 0
	
if args.nocoldate is not None:
	NoColDate = args.nocoldate
else:
	NoColDate = False
	
if args.nosubdate is not None:
	NoSubDate = args.nosubdate
else:
	NoSubDate = False

if args.nogisaid is not None:
	NoGisaid = args.nogisaid
else:
	NoGisaid = False
	
if args.noclade is not None:
	NoClade = args.noclade
else:
	NoClade = False

if args.noseqlen is not None:
	NoSeqLen = args.noseqlen
else:
	NoSeqLen = False
	
if args.maxloc is not None:
	MaxLocInfo = int(args.maxloc)
else:
	MaxLocInfo = 50
#	print("  Sequence length not defined. If you want to filter by sequence length, please use -seqlen option")	

### FUNCTIONS ###
# get date format with fields separated by "-"	
def getDateFmt (string):
	DateElements = string.split("-")
	if len(DateElements) == 3 and DateElements[0].isdigit() == True and DateElements[1].isdigit() == True and DateElements[2].isdigit()  == True:
		if len(DateElements[0]) == 2:
			Fmt = "%Y2d-%m-%d"
		elif len(DateElements[0]) == 4:
			Fmt = "%Y-%m-%d"
	elif len(DateElements) == 2 and DateElements[0].isdigit()  == True and DateElements[1].isdigit() == True:
		Fmt = "%Y-%m"
	elif DateElements[0].isdigit() == True:
		Fmt = "%Y"
	else:
		Fmt = "Invalid date format"

	return Fmt
		
# Get date frames
def GetDates (DateAfter, DateBefore):
	MyTimeDelta = datetime.timedelta(days = (DaysWindow)) 
	start = datetime.datetime.strptime(DateAfter, "%Y-%m-%d")
	curr = start
	if DateBefore == "":
		end = datetime.datetime.today()
	else:
		end = datetime.datetime.strptime(DateBefore, "%Y-%m-%d")
	NbrDays = (end - start).days
	TimeFrames = []
	for x in range (0, NbrDays, DaysWindow):
		if curr < (end - MyTimeDelta):
			TimeFrames.append([curr.strftime("%Y-%m-%d"), (curr + MyTimeDelta).strftime("%Y-%m-%d")])
			curr += MyTimeDelta + datetime.timedelta(days=1)
	
	TimeFramesDict = {}
	for framenbr in range (0, len(TimeFrames)): # insert time frame columns
		n = framenbr + 1
		FrameName = str(n) + "_" + TimeFrames[framenbr][0] + "_to_" + TimeFrames[framenbr][1]
		TimeFramesDict[FrameName] = [TimeFrames[framenbr], 0]  # for each FrameName: [[start_date, end_date], TotalGenomesInTheFrame-reserved space]
	
	return TimeFramesDict

# get comments from a file	
def getComments(FN):
	try:
		start = False
		my_splitlines = []
		InfoDict = {}
		extension = FN.split(".")[-1]
		if extension == "gz" or extension == "gzip":
			with gzip.open(FN, mode='rb') as f:
				for line in f:
					line = line.decode("utf-8")
					if line.split("|")[0][1:4] == 'END':	
						start = False
						break
					if start == True:
						my_splitlines.append(line)	
					if line.split("|")[0][1:6] == 'START':
						start = True
					
			f.close()
		else:
			with open(FN, mode='r') as f:
				for line in f:
					if line.split("|")[0][1:4] == 'END':	
						start = False
						break
					if start == True:
						my_splitlines.append(line)	
					if line.split("|")[0][1:6] == 'START':
						start = True
			f.close()

		for comment in my_splitlines:
			code = comment.split("|")[0][1:4]
			info = comment.split(":")[-1].strip()
			desc = comment.split("|")[-1].split(":")[0]
			InfoDict[code] = [desc, info]

		return InfoDict
		
	except FileNotFoundError:
                print("The file %s was not found !" % (FN))

# filter a table (pandas data frame)
def filterDF (DF):

        if MutCutOff != "":   
                DF = DF[(DF['#Mutations'] <= MutCutOff)]
                filters.append("NbrOfMut(< %i)" % MutCutOff)
                print("  Filtered genomes with less than %i mutations" % (MutCutOff))
        if DateAfter != "":   
                DF = DF[(DF['Collection_Date'] >= DateAfter)]
                filters.append("DateAfter(%s)" % DateAfter)
                print("  Filtered genomes with collection date after to %s" % (DateAfter))
        if DateBefore != "":   
                DF = DF[(DF['Collection_Date'] <= DateBefore)]
                filters.append("DateBefore(%s)" % DateBefore)
                print("  Filtered genomes with collection date before to %s" % (DateBefore))
        if Host != "":   
                DF = DF[(DF['Host'] == Host)]
                filters.append("Host(%s)" % Host)
                print("  Filtered genomes corresponding to %s host" % (Host))
        if Continent != "":   
                DF = DF[(DF['Continent'] == Continent)]
                filters.append("Continent(%s)" % Continent)
                print("  Filtered genomes corresponding to %s (continent)" % (Continent))
        if Country != "":   
                DF = DF[(DF['Country'] == Country)]
                filters.append("Country(%s)" % Country)
                print("  Filtered genomes corresponding to %s (country)" % (Country))
        if Region != "":   
                DF = DF[(DF['Region'] == Region)]
                filters.append("Region(%s)" % Region)
                print("  Filtered genomes corresponding to %s (region)" % (Region))
        if Clade != "":   
                DF = DF[(DF['Clade'] == Clade)]
                filters.append("Clade(%s)" % Clade)
                print("  Filtered genomes corresponding to %s (pango clade)" % (Clade))
        if SeqLen != "":   
                DF = DF[(DF['Sequence_length'] >= SeqLen)]
                filters.append("SeqLen(> %s bp)" % SeqLen)
                print("  Filtered genomes with sequence length >= %i bases" % (SeqLen))
        if Nbases != "":   
                DF = DF[(DF['#Ns'] <= Nbases)]
                filters.append("Nbases(%s)" % Nbases)
                print("  Filtered genomes with number of degenerated bases <= %i " % (Nbases))
        if MutList != "": 
        	for muta in MutList:
        		DF = DF[DF['PointMutations'].str.contains(muta)]
        		filters.append("PointMutations(%s)" % muta)
        		print("  Filtered genomes containing the mutations: %s" % (muta))
        		
        return DF

# filter a table (pandas data frame)
def end_filterDF (DF):

        if TopRank != 0:  
                if DF.shape[0] > TopRank:
                        DF = DF.iloc[0: TopRank - 1,]
                        print("  Filtered %i top ranked genomes" % (TopRank))
                        filters.append("Top ranked variants (< %i)" % TopRank)    
        if NoColDate != False:
                del DF["Collection_Date"]
                print("  Deleted collection date column")
                filters.append("No collection date output")
        if NoSubDate != False:
                del DF["Submission_Date"]
                print("  Deleted submission date column")
                filters.append("No submission date output")
        if NoGisaid != False:
                del DF["GISAID_ID"]
                print("  Deleted GISAID_ID column")
                filters.append("No GISAID ID output")
        if NoClade != False:
                del DF["Clade"]
                print("  Deleted Clade column")
                filters.append("No clade output")
        if NoSeqLen != False:
                del DF["Sequence_length"]
                print("  Deleted Sequence_length column")
                filters.append("No sequence length output")
 
        return DF


# get genotype related information
def getGenotypes(GenomeDF):
	## Get time frames
	DateBefore = InfoDict['GDD'][1]  # Assume the database date informed in the genome's table as comment as end date
	TimeFramesDict=GetDates(DateAfter, DateBefore) # Define time frames
	a = 0
	PercLast = 0
	TotalGenomes = GenomeDF.shape[0]
	GenotypesDict = {}
	for index, row in GenomeDF.iterrows():
		CollDate = row['Collection_Date']
		if row['PointMutations'] not in GenotypesDict.keys():
			GenotypesDict[row['PointMutations']] = {'Genotype' : row['PointMutations'],
					'GISAID_ID' : [row['GISAID_ID']],
					'#Mutations' : row['#Mutations'],
					'#Genomes' : 1,
					'Collection_Date' : [CollDate],
					'Submission_Date' : [row['Submission_Date']],
					'Continent' : {row['Continent'] : [1, 0]}, # [Dict key: values = counts, genotype frequency, TotalCounts]
					'Country' : {row['Country'] : [1, 0]}, # [Dict key: values = counts, genotype frequency, TotalCounts]
					'Region' : {row['Region'] : [1, 0]}, # [Dict key: values = counts, genotype frequency, TotalCounts]
					'Host' : [row['Host']],
					'Clade' : row['Clade'],
					'Sequence_length' : [row['Sequence_length']]
                        }
                        # assign 1 for the time frame corresponding to the collection date of a new genotype - 0 for all other time frames
			FormatOfCollDate = getDateFmt(CollDate)
			if FormatOfCollDate == "%Y-%m-%d" or FormatOfCollDate == "%Y2d-%m-%d": # if the collection date is in full format (2 or 4 digits year)
				if FormatOfCollDate == "%Y2d-%m-%d":
					DateElements = CollDate.split("-")
					DateElements[0] = '20' + DateElements[0]
					CollDate = "-".join(DateElements)
				FmtCollDate = datetime.datetime.strptime(CollDate, "%Y-%m-%d")
				for FrameName in TimeFramesDict.keys(): 
					StartDate = datetime.datetime.strptime(TimeFramesDict[FrameName][0][0], "%Y-%m-%d")
					EndDate = datetime.datetime.strptime(TimeFramesDict[FrameName][0][1], "%Y-%m-%d")
					if StartDate <= FmtCollDate and FmtCollDate <= EndDate: # insert time frame columns
						GenotypesDict[row['PointMutations']][FrameName] = 1
						TimeFramesDict[FrameName][1] += 1
					else:
						GenotypesDict[row['PointMutations']][FrameName] = 0	
			else: # if the collection date is not in full format do not assign
				for FrameName in TimeFramesDict.keys():
					GenotypesDict[row['PointMutations']][FrameName] = 0
                        
		else:  # update column values of the existing row
			GenotypesDict[row['PointMutations']]['GISAID_ID'].append(row['GISAID_ID'])
			GenotypesDict[row['PointMutations']]['#Genomes'] += 1
			GenotypesDict[row['PointMutations']]['Collection_Date'].append(row['Collection_Date'])
			GenotypesDict[row['PointMutations']]['Submission_Date'].append(row['Submission_Date'])
			if row['Continent'] not in GenotypesDict[row['PointMutations']]['Continent'].keys():
				GenotypesDict[row['PointMutations']]['Continent'][row['Continent']] = [1, 0]
			else:
				GenotypesDict[row['PointMutations']]['Continent'][row['Continent']][0] += 1			
			if row['Country'] not in GenotypesDict[row['PointMutations']]['Country'].keys():
				GenotypesDict[row['PointMutations']]['Country'][row['Country']] = [1, 0]
			else:
				GenotypesDict[row['PointMutations']]['Country'][row['Country']][0] += 1	
			if row['Region'] not in GenotypesDict[row['PointMutations']]['Region'].keys():
				GenotypesDict[row['PointMutations']]['Region'][row['Region']] = [1, 0]
			else:
				GenotypesDict[row['PointMutations']]['Region'][row['Region']][0] += 1	
			if row['Host'] not in GenotypesDict[row['PointMutations']]['Host']:
				GenotypesDict[row['PointMutations']]['Host'].append(row['Host'])
			if row['Sequence_length'] not in GenotypesDict[row['PointMutations']]['Sequence_length']:
				GenotypesDict[row['PointMutations']]['Sequence_length'].append(row['Sequence_length'])
			
			# add 1 to the time frame corresponding to the collection date of an existing genotype
			FormatOfCollDate = getDateFmt(CollDate)
			if FormatOfCollDate == "%Y-%m-%d" or FormatOfCollDate == "%Y2d-%m-%d": # if the collection date is in full format (2 or 4 digits year)
				if FormatOfCollDate == "%Y2d-%m-%d":
					DateElements = CollDate.split("-")
					DateElements[0] = '20' + DateElements[0]
					CollDate = "-".join(DateElements)
				FmtCollDate = datetime.datetime.strptime(CollDate, "%Y-%m-%d")
				for FrameName in TimeFramesDict.keys(): 
					StartDate = datetime.datetime.strptime(TimeFramesDict[FrameName][0][0], "%Y-%m-%d")
					EndDate = datetime.datetime.strptime(TimeFramesDict[FrameName][0][1], "%Y-%m-%d")
					if StartDate <= FmtCollDate and FmtCollDate <= EndDate:  # update time frame columns
						GenotypesDict[row['PointMutations']][FrameName] += 1
						TimeFramesDict[FrameName][1] += 1
				
		a += 1
		PercNow = int((a * 100)/TotalGenomes)
		if PercLast != PercNow:
			print("Creating genotype's table ... %i %s                                                    " % (PercNow, '%'), end="\r", flush=True)
			PercLast = PercNow

	print("Creating genotype's table ... Done!")
		
	print("Calculating location and time frame frequencies, slope and improving the output readability ... ", end="\r", flush=True)
	for genotype in GenotypesDict.keys():
		GenotypesDict[genotype]['GISAID_ID'] = ";".join(GenotypesDict[genotype]['GISAID_ID'])
		GenotypesDict[genotype]['Collection_Date'] = ";".join(GenotypesDict[genotype]['Collection_Date'])
		GenotypesDict[genotype]['Submission_Date'] = ";".join(GenotypesDict[genotype]['Submission_Date'])
		GenotypesDict[genotype]['Host'] = ";".join(GenotypesDict[genotype]['Host'])
		for n in range(len(GenotypesDict[genotype]['Sequence_length'])):
			GenotypesDict[genotype]['Sequence_length'][n] = str(GenotypesDict[genotype]['Sequence_length'][n])
		GenotypesDict[genotype]['Sequence_length'] = ";".join(GenotypesDict[genotype]['Sequence_length'])
		
		# calculate location related frequencies
		DataList = []
		GenotypesDict[genotype]['Continent'] = collections.OrderedDict(sorted(GenotypesDict[genotype]['Continent'].items(), key=lambda t: t[1][0], reverse=True)) # order by frequency	
		for item in GenotypesDict[genotype]['Continent'].keys():
			GenotypesDict[genotype]['Continent'][item][1] = GenotypesDict[genotype]['Continent'][item][0] / GenotypesDict[genotype]['#Genomes'] # calculate the genotype frequency in each cotinent
			if len(DataList) <= MaxLocInfo:
				if GenotypesDict[genotype]['Continent'][item][1] >= 0.01: # only report locations representing more than 1% of all genotypes (list length <= MaxLocInfo )
					data = item + ":" + str(GenotypesDict[genotype]['Continent'][item][0]) + ":" + str(GenotypesDict[genotype]['Continent'][item][1]) # more readable format
					DataList.append(data)
			else:
				break
		GenotypesDict[genotype]['Continent'] = ";".join(DataList)
		
		DataList = []
		GenotypesDict[genotype]['Country'] = collections.OrderedDict(sorted(GenotypesDict[genotype]['Country'].items(), key=lambda t: t[1][0], reverse=True)) # order by frequency
		for item in GenotypesDict[genotype]['Country'].keys():
			GenotypesDict[genotype]['Country'][item][1] = GenotypesDict[genotype]['Country'][item][0] / GenotypesDict[genotype]['#Genomes']  # calculate the genotype frequency in each country
			if len(DataList) <= MaxLocInfo:
				if GenotypesDict[genotype]['Country'][item][1] >= 0.01:  # only report locations representing more than 1% of all genotypes (list length <= MaxLocInfo )
					data = item + ":" + str(GenotypesDict[genotype]['Country'][item][0]) + ":" + str(GenotypesDict[genotype]['Country'][item][1]) # more readable format
					DataList.append(data)
			else:
				break
		GenotypesDict[genotype]['Country'] = ";".join(DataList)
		
		DataList = []
		GenotypesDict[genotype]['Region'] = collections.OrderedDict(sorted(GenotypesDict[genotype]['Region'].items(), key=lambda t: t[1][0], reverse=True)) # order by frequency
		for item in GenotypesDict[genotype]['Region'].keys():
			GenotypesDict[genotype]['Region'][item][1] = GenotypesDict[genotype]['Region'][item][0] / GenotypesDict[genotype]['#Genomes']  # calculate the genotype frequency in each country
			if len(DataList) <= MaxLocInfo:
				if GenotypesDict[genotype]['Region'][item][1] >= 0.01:  # only report locations representing more than 1% of all genotypes (list length <= MaxLocInfo )
					data = item + ":" + str(GenotypesDict[genotype]['Region'][item][0]) + ":" + str(GenotypesDict[genotype]['Region'][item][1]) # more readable format
					DataList.append(data)
			else:
				break
		GenotypesDict[genotype]['Region'] = ";".join(DataList)
		
		# calculate time frame frequencies and slope
		yall = []
		xall= []
		SlopeArray = []
		n = 1
		for FrameName in TimeFramesDict.keys():
			xall.append(n)  # create a list of values (number of values = number of collumns to be analyzed = number of time frames)
			n +=1
			if GenotypesDict[genotype][FrameName] > 0:
				GenotypesDict[genotype][FrameName] = GenotypesDict[genotype][FrameName] / TimeFramesDict[FrameName][1]
			yall.append(GenotypesDict[genotype][FrameName])
		y=[]
		for a in range (0, len(yall)-1): # take only values from the first increasing points to calculate slopes
			if yall[a+1] > yall[a]:
				y = yall[a+1:]
				x = xall[a+1:]
				break
		if len(y) > 1:
			pass
		else:
			y = yall
			x = xall

		if pd.isnull(linregress(x,y)[0]):
			GenotypesDict[genotype]['Slope'] = linregress([23,24], [yall[-2], yall[-1]])[0]
		else:
			GenotypesDict[genotype]['Slope'] = linregress(x, y)[0]
				
	print("Calculating location and time frame frequencies, slope and improving the output readability ... Done!")
	
	# Conversion of dictionary to DF
	GenotypeDF = pd.DataFrame.from_dict(GenotypesDict, orient='index')
	
	print("Ranking by slope ...                " , end="\r", flush=True)
	GenotypeDF = GenotypeDF.sort_values(by='Slope', ascending=False)
	print("Ranking by slope ... Done!               ")
	
	return GenotypeDF   

### MAIN WORK ###
## Verbose section
# Starting screen follow-up
print("\nRunning ...")
LogFile=open("Genotypes.log", "w") # for logging in future implementations

# Read genomes table and comments
print("\nReading input genome's table (file: %s) ...        " % (GenTableFN), end="\r", flush=True)
try:
	InfoDict = getComments(GenTableFN)
	skiplines = len(InfoDict.keys()) + 2
	extension = GenTableFN.split(".")[-1]
	if  extension == "gz" or extension == "gzip":
		GenomeDF = pd.read_csv(GenTableFN, sep='\t', compression='gzip', skiprows=skiplines)
	else:
		GenomeDF = pd.read_csv(GenTableFN, sep='\t', skiprows=skiplines)
	
	print("\nReading input genome's table (file: %s) ... Done!       " % (GenTableFN))
except FileNotFoundError:
	print("The file %s was not found !" % (filename))

# Filter following user provided options
print("Filtering by default values and/or user options:")
UnfilteredNbrOfGeno = GenomeDF.shape[0]
print(" Number of unfiltred genomes: %i" % (UnfilteredNbrOfGeno))
GenomeDF = filterDF(GenomeDF) 
FilteredNbrOfGeno = GenomeDF.shape[0]  
print(" Number of filtred genomes: %i" % (FilteredNbrOfGeno))

# Gather genotype data using the filtered data frame
GenotypeDF = getGenotypes(GenomeDF)
del GenomeDF

# Final filtering
print("Applying final filter ...                ")
EndUnfilteredNbrOfGeno = GenotypeDF.shape[0]
GenotypeDF = end_filterDF(GenotypeDF) 
EndFilteredNbrOfGeno = GenotypeDF.shape[0]  
print("  Final number genomes reported: %i" % (EndFilteredNbrOfGeno))
print("Applying final filter ... Done!               ")

# Insert the rank column	
Rank = [i for i in range(1, GenotypeDF.shape[0] + 1)]
GenotypeDF.insert(loc = 0, column = 'Rank', value = Rank )
del Rank

print("Sample of the results - currently top emerging genotypes:\n")
print(GenotypeDF)

# Write table
print("Writing genotype's table: %s" % (OutFN), end="\r", flush=True)
with gzip.open(OutFN, 'wt') as f:
	f.write("#START| **** INFORMATION SECTION ****\n")
	for code in InfoDict.keys():
		f.write("#%s|%s: %s\n" % (code, InfoDict[code][0], InfoDict[code][1]))
	f.write("#NFG| Number of filtered genomes: %s\n#AUF| Used filters: %s\n#TFD| Time frame(days): %i\n#END| **** INFORMATION SECTION ****\n" % (FilteredNbrOfGeno, ",".join(filters), DaysWindowDef))
	GenotypeDF.to_csv(f, index = False, sep="\t", mode='a')
print("Writing genotype's table: %s, done!" % (OutFN))

LogFile.close()

print("\nJob finished !" )
quit()

