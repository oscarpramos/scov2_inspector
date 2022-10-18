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
GenTableFN = OutFN = MutCutOff = MutaList = Host = DateBefore = DateAfter = DaysWindow = Continent = Country = Region = Clade = SeqLen = Nbases = ""
filters = []
    
### INITIALIZING TASKS ###
# Help about how to use the script
parser = argparse.ArgumentParser(description="01b_getSNPTable.py: Create a table containing snp related data using Genome's table as input. \n Usage: ./01b_getSNPTable..py -i <GenomesTable(.tsv or .tsv.gz)> <option>\n", prefix_chars='-+')
parser.add_argument("-i", help="File name of the table containing all mutatated genomes and related information (.tsv or tsv.gz)")
parser.add_argument("-o", help="Output file name (tab separated)")
parser.add_argument("-mutN", help="Filer based on maximum number of mutations per genotype  (default = 100)")
parser.add_argument("-mutL", help="Filter mutations using a list of mutations separated by ',' (default = none; nucleotide, protein and annotations allowed")
parser.add_argument("-host", help="Filter by host name (default = Human)")
parser.add_argument("-da", help="Filter by collection date after (YYYY-MM-DD) (default = none)")
parser.add_argument("-db", help="Filter by collection date before (YYYY-MM-DD) (default = none)")
parser.add_argument("-cont", help="Filter by continent name (default = none)")
parser.add_argument("-country", help="Filter by country name (default = none)")
parser.add_argument("-region", help="Filter by region name (default = none)")
parser.add_argument("-clade", help="Filter by pango clade (default = none)")
parser.add_argument("-seqlen", help="Filter by minimum sequence length (default = 29000)")
parser.add_argument("-nbases", help="Filter by the maximum number of degenerated bases found in sequence (default = 704)")

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
        OutFN = "IndividualGenomeList_filtered.tsv.gz"

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

### FUNCTIONS ###
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
                filters.append("NbrOfMut(%i)" % MutCutOff)
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
                filters.append("SeqLen(%s)" % SeqLen)
                print("  Filtered genomes with sequence length >= %i " % (SeqLen))
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
        
### MAIN WORK ###
## Verbose section
# Starting screen follow-up
print("\nRunning ...")
LogFile=open("SNPs.log", "w") # for logging in future implementations

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
	
	print("Reading input genome's table (file: %s) ... Done!       " % (GenTableFN))
except FileNotFoundError:
	print("The file %s was not found !" % (filename))

# Filter following user provided options
print("Filtering by default values and/or user options:")
UnfilteredNbrOfGeno = GenomeDF.shape[0]
print(" Number of unfiltred genomes: %i" % (UnfilteredNbrOfGeno))
GenomeDF = filterDF(GenomeDF) 
FilteredNbrOfGeno = GenomeDF.shape[0]  
print(" Number of filtred genomes: %i" % (FilteredNbrOfGeno))

# Gather snp data using the filtered data frame
GenomeDF['Clade'] = GenomeDF['Clade'].replace(np.nan, "Not known")

# Write table
print("Writing genotype's table: %s" % (OutFN), end="\r", flush=True)
with gzip.open(OutFN, 'wt') as f:
	f.write("#START| **** INFORMATION SECTION ****\n")
	for code in InfoDict.keys():
		f.write("#%s|%s: %s\n" % (code, InfoDict[code][0], InfoDict[code][1]))
	f.write("#NFG| Number of filtered genomes: %s\n#AUF| Used filters: %s\n#END| **** INFORMATION SECTION ****\n" % (FilteredNbrOfGeno, ",".join(filters)))
	GenomeDF.to_csv(f, index = False, sep="\t", mode='a')
print("Writing snp's table: %s, done!" % (OutFN))
LogFile.close()
