#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### INFO ###
__author__ = "Oscar Ramos"
__copyright__ = "Copyright 20120, Oscar Ramos"
__license__ = ""
__maintainer__ = "Oscar Ramos"
__version__ = "0.1"
__email__ = "oscar.pereira-ramos@cea.fr"
__status__ = "alpha"

### MODULES IMPORT ###
import argparse
import gzip
import collections
import locale
import warnings
import pandas as pd
import primer3
import plotly.express as px
import plotly.graph_objects as go

## American global settings
locale.setlocale(locale.LC_ALL, 'C')

### GLOBAL VARIABLES DEFINITION ###
GenTableFN = OutFN = MutCutOff = MutaList = Host = DateBefore = DateAfter = DaysWindow = Continent = Country = Region = Clade = SeqLen = ""
filters = []
    
### INITIALIZING TASKS ###
# Help about how to use the script
parser = argparse.ArgumentParser(description="02b_getSNPTable.py: Create a table of conservation by position and evaluate qPCR sets using Genome's table as input. \n Usage: ./02c_getConsTable+qPCR_2gz..py -i <SNPTable(.tsv or .tsv.gz)> <options>\n", prefix_chars='-+')
parser.add_argument("-i", help="File name of the table containing all mutatated genomes and related information (.tsv or tsv.gz)")
parser.add_argument("-q", help="The tab separated file containing qPCR primers and probe information")
parser.add_argument("-r", help="The reference genome (fasta file)")
parser.add_argument("-g", help="The gff file containing the annotation for the reference genome (gff)")
parser.add_argument("-o", help="Output file name (tab separated)")
parser.add_argument("-mutN", help="Filer based on maximum number of mutations per genotype  (default = 100)")
parser.add_argument("-mutL", help="Filter mutations using a list of mutations separated by ',' (default = none; nucleotide, protein and annotations allowed)")
parser.add_argument("-host", help="Filter by host name (default = Human)")
parser.add_argument("-da", help="Filter by collection date after (YYYY-MM-DD) (default = none)")
parser.add_argument("-db", help="Filter by collection date before (YYYY-MM-DD) (default = none)")
parser.add_argument("-cont", help="Filter by continent name (default = none)")
parser.add_argument("-country", help="Filter by country name (default = none)")
parser.add_argument("-region", help="Filter by region name (default = none)")
parser.add_argument("-clade", help="Filter by pango clade (default = none)")
parser.add_argument("-seqlen", help="Filter by minimum sequence length (default = 29000)")
parser.add_argument("-nbases", help="Filter by the maximum number of degenerated bases found in sequence (default = 704)")
parser.add_argument("-t", help="Temperature in °C (default: 37)")
parser.add_argument("-mv", help="monovalent ions concentration in mM (default: 50)")
parser.add_argument("-dv", help="divalent ions concentration in mM (default: 1.5)")
parser.add_argument("-dntp", help="dNTP concentration in mM (default: 0.3)")
parser.add_argument("-dna", help="DNA oligonucleotide concentration in nM (default: 500)")

#actually parse command line options 
args = parser.parse_args()

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
InfoDict = getComments(args.i)

if args.i is not None:
	GenTableFN = args.i
else:
	print("  Please provide genomes in fasta format provided by GISAID (-i)")
	quit()
	
if args.q is not None:
	qPCRFN = args.q
else:
	print("  Please provide qPCR primers/probes information in .tsv format (-q)")
	quit()
	
if args.r is not None:
	RefGenoFN = args.r
else:
	try:
		RefGenoFN = InfoDict['RFN'][1]
	except:
		print("  Please provide reference genome in fasta format (-r)")
		quit()
	
if args.g is not None:
	GffFN = args.g
else:
	try:
		GffFN = InfoDict['AFN'][1]
	except:
		print("  Please provide gff annotation file (-g)")
		quit()

if args.o is not None:
	OutFN = args.o
	extension = args.o.rsplit(".")[-1]
	if extension != "gz" and  extension != "gzip":
	      OutFN = OutFN + '.gz'  
else:
        OutFN = "ConservationTable.tsv"

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

if args.t is not None:
        Temp = int(args.t)
else:
        Temp = 37

if args.dntp is not None:
        dNTPConc = int(args.dntp)
else:
	dNTPConc = 0.3

if args.mv is not None:
        MVConc = args.mv
else:
        MVConc = 60

if args.dv is not None:
        DVConc = args.dv
else:
        DVConc = 1.5

if args.dna is not None:
        DNAConc = args.dna
else:
        DNAConc = 500

### FUNCTIONS ###
# transfer the content of a file to the RAM  
def get_file_content(filename):
	try:
		with open(filename) as f:
		        data = f.read()
		        my_splitlines = data.splitlines()
		        f.close()
		        return my_splitlines
	except FileNotFoundError:
		print("The file %s was not found !" % (filename))
                
# transfer the content of a file to the RAM  
def get_gzfile_content(filename):
	try:
		with gzip.open(filename, mode='rb') as f:
			data = f.read().decode("utf-8")
			my_splitlines = data.splitlines()
			f.close()
			return my_splitlines
	except FileNotFoundError:
		print("The file %s was not found !" % (filename))  
		

# parse qPCR information file (.tsv)              
def parse_qPCRdata(file):
	qPCRdata = []
	qPCR_tsv = get_file_content(file)
	for line in qPCR_tsv:
                        Fields = line.rsplit("\t")
                        qPCRdata.append(Fields)
	return qPCRdata
              
# load FastaDict from fasta file	
def getFastaData (FastaFileName):
	extension = FastaFileName.split(".")[-1]
	if  extension == "gz" or extension == "gzip":
		fasta = get_gzfile_content(FastaFileName)
	else:
		fasta = get_file_content(FastaFileName)
	a = 0
	TotalLines = len(fasta)
	PercLast = 0
	header = ""
	seq = ""
	FastaDict = {}
	FirstSeqLine = True
	for line in fasta:
		if line[:1] == ">":
			header = line.replace(">", "").strip()
			FirstSeqLine = True
			a += 1
		else:
			if FirstSeqLine == True:
				seq = line.upper().replace("\n","")
				FirstSeqLine = False
			else:
				seq = "%s\n%s" % (seq, line.upper().replace("\n","")) 
		FastaDict[header] = seq
		a += 1
		PercNow = int((a * 100)/TotalLines)
		if PercLast != PercNow:
			print("Reading genomes from file %s ...                    " % (FastaFileName), end="\r", flush=True)
			PercLast = PercNow
	print("Reading genomes from file %s ... Done!                    " % (FastaFileName))
		
	return a, FastaDict

# Get regions from a GFF file 
def GFF2CDSs(GFFdata):
        CDSs = []
        try:
                GFFdata = get_file_content(GffFN) # Read GFF file
                for line in GFFdata:
                        Fields = line.rsplit("\t")
                        if Fields[2] == "Mature_peptide" or Fields[2] == "UTR": # Get mature peptide and UTR annotations
                                CDSs.append([Fields[3], Fields[4], Fields[8].rsplit(";")[0]])  # Start, Stop, CDS name         
   
                return CDSs

        except ValueError:
                return "Could not get CDSs!"

# Verify if a mutation occurs inside a region (region extracted from a provided GFF file) 
def IFCDS (CDSs, now):
	now = int(now)
	CDSHit = []
	for CDS in CDSs:
		if int(CDS[0]) <= now <= int(CDS[1]):
			if CDSHit == []:
				CDSHit = [CDS[2].strip()]
			else:
				CDSHit.append(CDS[2].strip())
			
	CDSHit = "_&_".join(CDSHit)
	    	
	return CDSHit
        
# calculate percentage of G+C 
def perc_gc(seq):
	x = 0
	seq = seq.upper()
	for char in seq:
		if char == "G" or char == "C":
			x += 1
	PercGC = (x *100) / len(seq)
	return PercGC
	
# get the reverse complement of a DNA sequence
def rc_dna(sequence):
	
	rc_sequence = ''
	reverse = sequence[::-1]
	
	# complementary bases dictionary
	complement_nt = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N' }

	for n in range(0, len(reverse)):
		if reverse[n] in complement_nt.keys():
			rc_sequence += complement_nt.get(reverse[n], "N")

	return rc_sequence
	
# retrieve the minimlal length of a array of oligonucleotides
def getMinLen(qPCRData, Colnbr):
	MinLen = 10000
	for SetNbr in range(1, len(qPCRData)):
		if len(qPCRData[SetNbr][Colnbr]) < MinLen:
			MinLen = len(qPCRData[SetNbr][Colnbr])
	return MinLen

# retrieve the minimlal length of a array of oligonucleotides
def getMaxLen(qPCRData, Colnbr):
	MaxLen = 0
	for SetNbr in range(1, len(qPCRData)):
		if len(qPCRData[SetNbr][Colnbr]) > MaxLen:
			MaxLen = len(qPCRData[SetNbr][Colnbr])
	return MaxLen

# Translate the consequence of NT mutations to aa (Special CDS with ribosome slippage; i.e., RdRP CDS in the case of SARS-CoV-2)
def get_RdRP_Codonmut(Muta, CDS):  # Normal translation until codon 9, Ribosome splippage Nt 27
	lag = 1
	LastNormalCodon = 9
	MutaPos = int(Muta[1:-1])
	RefNT = Muta[0]
	QueryNT = Muta[-1]
	CDS_start = int(CDS[0]) - 1
	CDS_stop = int(CDS[1])
	CDSSeq = RefGenomeSeq[CDS_start:CDS_stop]
	if ((MutaPos-CDS_start) / 3)%1 == 0:
		codon_info = [int((((MutaPos)-CDS_start) / 3 )) -1 , 2] # codon and position
	if 0.4 < ((MutaPos-CDS_start) / 3)%1 < 0.8:
		codon_info = [int(((MutaPos)-CDS_start) / 3), 1]
	if 0 < ((MutaPos-CDS_start) / 3)%1 < 0.4:
		codon_info = [int(((MutaPos)-CDS_start) / 3), 0]
	for n in range(0,len(CDSSeq),3):
		codon_nbr = int(((n+3) / 3) -1)
		if codon_nbr < LastNormalCodon:
			if codon_info[0] == (codon_nbr):
				old_codon = CDSSeq[n:n+3]
				new_codon = list(old_codon)
				new_codon[codon_info[1]] = QueryNT
				new_codon = str(''.join(new_codon))
				old_aa = translate_dna(old_codon)
				new_aa = translate_dna(new_codon)
				if old_aa != new_aa:
					AAresult = old_aa + str(codon_nbr + 1) + new_aa
				else:
					AAresult = old_aa + str(codon_nbr + 1) + old_aa
		else:
			n = n - lag
			if codon_info[1] != 2:
				if codon_info[0] == (codon_nbr):
					old_codon = CDSSeq[n:n+3]
					new_codon = list(old_codon)
					new_codon[codon_info[1] + lag] = QueryNT
					new_codon = str(''.join(new_codon))
					old_aa = translate_dna(old_codon)
					new_aa = translate_dna(new_codon)
					if old_aa != new_aa:
						AAresult = old_aa + str(codon_nbr + 1) + new_aa
					else:
                                		AAresult = old_aa + str(codon_nbr + 1) + old_aa

			else:
				if (codon_info[0] + 1) == codon_nbr :
					old_codon = CDSSeq[n:n+3]
					new_codon = list(old_codon)
					new_codon[0] = QueryNT
					new_codon = str(''.join(new_codon))
					old_aa = translate_dna(old_codon)
					new_aa = translate_dna(new_codon)
					if old_aa != new_aa:
						AAresult = old_aa + str(codon_nbr + 1) + new_aa
					else:
						AAresult = old_aa + str(codon_nbr + 1) + old_aa

	return AAresult
	del CDS   

# Get codon info from RDRP CDS
def get_RDRPCodonInfo(RefGenomeSeq, Posnbr, CDS_start, CDS_stop):
	codonnbr = ((Posnbr - CDS_start) / 3) + 1
	SlippagePos = 13468
	SlippageCodonPos = 3
	SlippageLag = 1
	codon_info = []
	if Posnbr < SlippagePos: # expected reading frame region
		codonnbr = ((Posnbr - CDS_start) / 3) + 1
		CodonPosFactor = ((Posnbr + 1) - CDS_start ) / 3
		if CodonPosFactor %1 == 0: 
			codon_info = [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr, 3] # codonnbr and positionnbr
		if 0.4 < CodonPosFactor %1 < 0.8:
			codon_info = [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
		if 0 < CodonPosFactor %1 < 0.4:
			codon_info = [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr, 1]
	elif  Posnbr == SlippagePos: # 2 codons are changed at this position
		codonnbr = ((Posnbr - CDS_start) / 3) + 1
		CodonPosFactor = ((Posnbr + 1) - CDS_start ) / 3
		if CodonPosFactor %1 == 0:
			codon_info = [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr, 3] # codonnbr and positionnbr
		if 0.4 < CodonPosFactor %1 < 0.8:
			codon_info = [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
		if 0 < CodonPosFactor %1 < 0.4:
			codon_info = [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr, 1]	
	else: # shifted reding frame
		codonnbr = ((Posnbr - CDS_start) / 3) + 1
		CodonPosFactor = ((Posnbr + 1) - CDS_start ) / 3
		if CodonPosFactor %1 == 0:
			codon_info = [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr, 3] # codonnbr and positionnbr
		if 0.4 < CodonPosFactor %1 < 0.8:
			codon_info = [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
		if 0 < CodonPosFactor %1 < 0.4:
			codon_info = [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr, 1]
		
	
#	input('Posnbr: %i\tCodon: %s\tCodonnbr: %i\tPosInCodon: %i' % (Posnbr, codon_info[0], codon_info[1], codon_info[2]))
	return codon_info

# Get codon info from standard CDS
def get_CodonInfo(RefGenomeSeq, Posnbr, CDS_start, CDS_stop):
	codon_info = []
	codonnbr = ((Posnbr - CDS_start) / 3) + 1
	CodonPosFactor = ((Posnbr + 1) - CDS_start ) / 3
	if CodonPosFactor %1 == 0:
		codon_info = [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr, 3] # codonnbr and positionnbr
	if 0.4 < CodonPosFactor %1 < 0.8:
		codon_info = [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
	if 0 < CodonPosFactor %1 < 0.4:
		codon_info = [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr, 1]
	
#	input('Posnbr: %i\tCodon: %s\tCodonnbr: %i\tPosInCodon: %i' % (Posnbr, codon_info[0], codon_info[1], codon_info[2]))
	return codon_info

def getSilentProbab (RefGenomeSeq, ):

	# dictionary to retrieve codons for each AA
	AA1Lett2CodonDict = { 'A' : ('GCT', 'GCC', 'GCA', 'GCG'),
	'B' : ('GAT', 'GAC', 'AAT', 'AAC'), # D or N
	'C' : ('TGT', 'TGC'),
	'D' : ('GAT', 'GAC'),
	'E' : ('GAA', 'GAG'),
	'F' : ('TTT', 'TTC'),
	'G' : ('GGT', 'GGC', 'GGA', 'GGG'),
	'H' : ('CAT', 'CAC'),
	'I' : ('ATT', 'ATC', 'ATA'),
	'K' : ('AAA', 'AAG'),
	'L' : ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
	'M' : ('ATG'),
	'N' : ('AAT', 'AAC'),
	'P' : ('CCT', 'CCC', 'CCA', 'CCG'),
	'Q' : ('CAA', 'CAG'),
	'R' : ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
	'S' : ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
	'T' : ('ACT', 'ACC', 'ACA', 'ACG'),
	'V' : ('GTT', 'GTC', 'GTA', 'GTG'),
	'W' : ('TGG'),
	'X' : ('TAA', 'TAG', 'TGA'),
	'Y' : ('TAT', 'TAC'),
	'Z' : ('GAA', 'GAG', 'CAA', 'CAG') # E or Q
	}

# return a position vs conservation table	
def getConservArray (RefGenomeSeq, SNPsDict, FilteredNbrOfGeno, CDSs):
	CDSsDict = { CDS[2].strip() :  (CDS[0], CDS[1]) for CDS in CDSs }
	ConservArray = [None] * len(RefGenomeSeq)
	PosArray = [None] * len(RefGenomeSeq)
	RefNTArray = [None] * len(RefGenomeSeq)
	AnnotArray = [None] * len(RefGenomeSeq)
	NtmutArray = [None] * len(RefGenomeSeq)
	AAmutArray = [None] * len(RefGenomeSeq)
	NTAccRateArray = [None] * len(RefGenomeSeq)
        

	for i in range(0, len(RefGenomeSeq)):
		NTList = []
		AAList = []
		NTAccRateN = 0
		NTAccRateS = 0
		AccRate = 0
		Posnbr = i + 1
		PosArray[i] = Posnbr
		if Posnbr in SNPsDict.keys():
			NonConservSum = 0
			for NTQuery in SNPsDict[Posnbr]['NTQuery']:
				NonConservSum += SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				NTList.append(SNPsDict[Posnbr]['NTRef'] + str(Posnbr) + NTQuery + "(" + "{:f}".format(SNPsDict[Posnbr]['NTQuery'][NTQuery][2]) + ")")
				AAList.append(SNPsDict[Posnbr]['NTQuery'][NTQuery][1] + "(" + "{:f}".format(SNPsDict[Posnbr]['NTQuery'][NTQuery][2]) + ")")
				SkipUTR = ("5'UTR", "3'UTR","UTR")
				if SNPsDict[Posnbr]['Annotation'] not in SkipUTR:
					if SNPsDict[Posnbr]['NTQuery'][NTQuery][1][0] == SNPsDict[Posnbr]['NTQuery'][NTQuery][1][-1]:
						NTAccRateS +=  SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
					else:
						NTAccRateN += SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
					if SNPsDict[Posnbr]['Annotation'] != 'nsp12(RdRP)':
#						print(CDSsDict)
#					print('Posnbr = %i\tNTRef = %s' % (Posnbr, SNPsDict[Posnbr]['NTRef']) )
						codon_info = get_CodonInfo(RefGenomeSeq, Posnbr, int(CDSsDict[SNPsDict[Posnbr]['Annotation']][0]), int(CDSsDict[SNPsDict[Posnbr]['Annotation']][1]))
					else:
						codon_info = get_RDRPCodonInfo(RefGenomeSeq, Posnbr, int(CDSsDict[SNPsDict[Posnbr]['Annotation']][0]), int(CDSsDict[SNPsDict[Posnbr]['Annotation']][1]))
				else:
					AccRate = '-'
			
			NtmutArray[i] = ";".join(NTList)
			AAmutArray[i] = ";".join(AAList)
			RefNTArray[i] = SNPsDict[Posnbr]['NTRef']
			ConservArray[i] = (1 - NonConservSum)
			AnnotArray[i] = SNPsDict[Posnbr]['Annotation']
			# Calculate AccRatio
			
			if NTAccRateS == NTAccRateN == 0 :
				NTAccRatio = 1.0
			elif NTAccRateS == 0:
				NTAccRatio = '>1'
			elif NTAccRateN == 0:
				NTAccRatio = '<1'
			else:
				NTAccRatio = NTAccRateN / NTAccRateS
			NTAccRateArray[i] = NTAccRatio

			
		else:
			NtmutArray[i] = "fully conserved"
			AAmutArray[i] = "fully conserved"
			RefNTArray[i] = RefGenomeSeq[i]
			NTAccRateArray[i] = "fully conserved"
			ConservArray[i] = 1
			AnnotArray[i] = IFCDS(CDSs, Posnbr) # get annotation associated to the mutation
			
		
		
	df = pd.DataFrame({'PosNbr' : PosArray, 'RefNT' : RefNTArray, 'Conservation' : ConservArray, 'NTmutation(s)' : NtmutArray, 'AAmutation(s)' : AAmutArray, 'Annotation(s)' : AnnotArray, 'AcceptanceRate (w)' : NTAccRateArray})
	
	return df

# calculate frequencies of mutations for each position 
def getFinalPosFreqs (SNPsDict, Posnbr, strand):
	Afreq = Cfreq = Gfreq = Tfreq = Delfreq = InsFreq = NFreq = 0.00
	if strand != "-":
		if Posnbr in SNPsDict.keys():
			for NTQuery in SNPsDict[Posnbr]['NTQuery'].keys():
				if NTQuery == "A":
					Afreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == "C":
					Cfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == "G":
					Gfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == "T":
					Tfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == ".":
					Delfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				else:
					NFreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
			if SNPsDict[Posnbr]['NTRef'] ==  "A":
				Afreq = 1 - (Cfreq + Gfreq + Tfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  "C":
				Cfreq = 1 - (Afreq + Gfreq + Tfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  "G":
				Gfreq = 1 - (Afreq + Cfreq + Tfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  "T":
				Tfreq = 1 - (Afreq + Cfreq + Gfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  ".":
				Insfreq = 1 - (Afreq + Cfreq + Gfreq + Tfreq + NFreq)
		else:
			if RefGenomeSeq[Posnbr -1 : Posnbr] ==  "A":
				Afreq = 1
			elif RefGenomeSeq[Posnbr -1 : Posnbr] ==  "C":
				Cfreq = 1
			elif RefGenomeSeq[Posnbr -1 : Posnbr] ==  "G":
				Gfreq = 1
			elif RefGenomeSeq[Posnbr -1 : Posnbr] ==  "T":
				Tfreq = 1
			Delfreq = 0
			InsFreq = 0
	else:
		if Posnbr in SNPsDict.keys():
			for NTQuery in SNPsDict[Posnbr]['NTQuery'].keys():
				if NTQuery == "A":
					Tfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == "C":
					Gfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == "G":
					Cfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == "T":
					Afreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				elif NTQuery == ".":
					Delfreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
				else:
					NFreq = SNPsDict[Posnbr]['NTQuery'][NTQuery][2]
			if SNPsDict[Posnbr]['NTRef'] ==  "A":
				Tfreq = 1 - (Cfreq + Gfreq + Tfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  "C":
				Gfreq = 1 - (Afreq + Gfreq + Tfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  "G":
				Cfreq = 1 - (Afreq + Cfreq + Tfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  "T":
				Afreq = 1 - (Afreq + Cfreq + Gfreq + Delfreq + NFreq)
			elif SNPsDict[Posnbr]['NTRef'] ==  ".":
				Insfreq = 1 - (Afreq + Cfreq + Gfreq + Tfreq + NFreq)
		else:
			if RefGenomeSeq[Posnbr -1 : Posnbr] ==  "A":
				Tfreq = 1
			elif RefGenomeSeq[Posnbr -1 : Posnbr] ==  "C":
				Gfreq = 1
			elif RefGenomeSeq[Posnbr -1 : Posnbr] ==  "G":
				Cfreq = 1
			elif RefGenomeSeq[Posnbr -1 : Posnbr] ==  "T":
				Afreq = 1
			Delfreq = 0
			InsFreq = 0

	return Afreq, Cfreq, Gfreq, Tfreq, Delfreq, InsFreq, NFreq

# get comments from a file	
def getComments(FN):
	try:
		start = False
		my_splitlines = []
		InfoDict = {}
		extension = GenTableFN.split(".")[-1]
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

# get snp related information
def getSNPs(GenomeDF, RefGenomeSeq):
	a = 0
	PercLast = 0
	TotalGenomes = GenomeDF.shape[0]
	SNPsDict = {}
	for index, row in GenomeDF.iterrows():
		MutasList = row['PointMutations'].rsplit("|")
		for muta in MutasList:
				NT = muta.split(";")[0]
				NTPos = int(NT[1:-1])
				NTRef = RefGenomeSeq[NTPos - 1]
				NTQuery = NT[-1]
				AA = muta.split(";")[1]
				Annot = muta.split(";")[2]
				# create a dictionary using positions as keys
				if NTPos not in SNPsDict.keys():
					SNPsDict[NTPos] = {'NTRef' : NTRef,
						'NTQuery' : {NTQuery : [1, AA, 0]}, # number of genomes, aa mutation, frequency of the mutation (initialization using 0) 
						'Annotation' : Annot 
		                		}
				else:
					if NTQuery not in SNPsDict[NTPos]['NTQuery'].keys():
						SNPsDict[NTPos]['NTQuery'][NTQuery] = [1, AA, 0]  # number of genomes, aa mutation, annotation and frequency of the mutation (initialization using 0) 
					else:
						SNPsDict[NTPos]['NTQuery'][NTQuery][0] += 1
		a += 1
		PercNow = int((a * 100)/TotalGenomes)
		if PercLast != PercNow:
			print("Creating snp's table ... %i %s                                                    " % (PercNow, '%'), end="\r", flush=True)
			PercLast = PercNow
	
	print("Creating snp's table ... Done!")
		
	a = 0
	for NTPos in SNPsDict.keys():
		if PercLast != PercNow:
			PercLast = PercNow
		for NTChange in SNPsDict[NTPos]['NTQuery'].keys():
			SNPsDict[NTPos]['NTQuery'][NTChange][2] = SNPsDict[NTPos]['NTQuery'][NTChange][0] / TotalGenomes
		a += 1
		print("Calculating mutation frequencies ... %i %s                                                    " % (PercNow, '%'), end="\r", flush=True)
		PercNow = int((a * 100)/TotalGenomes)
	print("Calculating mutation frequencies ... Done!                                                    ", end="\r", flush=True) 
                
	SNPsDict = collections.OrderedDict(sorted(SNPsDict.items(), key=lambda t: t[0], reverse=False)) # only to order based on genome position
	for Posnbr in SNPsDict.keys():
		SNPsDict[Posnbr]['NTQuery'] = collections.OrderedDict(sorted(SNPsDict[Posnbr]['NTQuery'].items(), key=lambda t: t[1][2], reverse=True)) # order mutations for each position based on their frequencies
	
	return SNPsDict

# Analyzed qPCR oligos using information of conservation for each position in reference genome
def do_qPCR_analysis(qPCRData, ConservDF, RefGenomeSeq):             
        qPCResults = {}
        for qPCRset in range(1,len(qPCRData)):
	        SetNbr = qPCRData[qPCRset][0]   
	        Name = "%s_%s_%s_%s" % (qPCRData[qPCRset][0], qPCRData[qPCRset][1], qPCRData[qPCRset][2], qPCRData[qPCRset][3]) 
	        REM = "%s; %s; %s" % (qPCRData[qPCRset][1], qPCRData[qPCRset][2], qPCRData[qPCRset][3])
	        qPCRsetOutFile=open((Name + ".out"), "w")
	        qPCRsetOutFile.write("NAME qPCR_set_%s\n# REMARKS\nREM %s\n#Number of analyzed genomes\nNbrG %i\n" % (qPCRData[qPCRset][0], REM, FilteredNbrOfGeno))
	        qPCResults[SetNbr] = [Name, FilteredNbrOfGeno, 0, {}] # qPCRSetNumber : Name, NoG, TotalScore, OligosDict
	        # cycle over all oligonucleotides in the qPCRset 	
	        for oligoname in oligotype: 
		        # Wipe spaces in oligonucleotide sequences
		        Seq = qPCRData[qPCRset][oligotype[oligoname][0]].strip()
		        Consider = False
		        primerindex = 1
		        Tm = primer3.calcTm(Seq, mv_conc = MVConc, dv_conc = DVConc, dntp_conc=dNTPConc, dna_conc=DNAConc, tm_method='santalucia', salt_corrections_method="owczarzy")
		        dG = (primer3.bindings.calcHeterodimer(Seq, rc_dna(Seq), mv_conc = MVConc, dv_conc = DVConc, dntp_conc=dNTPConc, dna_conc=DNAConc, temp_c=Temp, max_loop=30, output_structure=False)).dg
		        GC = perc_gc(Seq)
		        LN = len(Seq)
		        qPCRsetOutFile.write("# Section: %s oligo \nEvery line starting by %s is for %s oligo\n# %sGC : %s of GC \n%sGC %0.2f\n# %sLN : length in bases\n%sLN %s\n# %sTM : Tm (degrees Celcius)\n%sTM %0.2f\n# %sDG: Free energy (dG) for perfect/complete complementarity (kcal/mol)\n%sDG %0.2f\n" % (oligoname, oligotype[oligoname][3], oligoname, oligotype[oligoname][3], '%', oligotype[oligoname][3], GC, oligotype[oligoname][3], oligotype[oligoname][3], LN, oligotype[oligoname][3], oligotype[oligoname][3], Tm, oligotype[oligoname][3], oligotype[oligoname][3], dG))
		        qPCRsetOutFile.write("# %sN tag the lines with nucleotide statistics\n# col1=tag_%sN+relative_index;\n# col2=nucleotide_in_primer;\n# col3=absolute_index_position;\n# col4-7=Conservation frequency of A, C, G, T, respectively, present in the reference genome. \n# col8-9=Frequency of deletions and insertion, respectively, observed in the compared genomes\n# colx=Other features ?...\n# Oligo_info\tRefPos\tConservA\tConservC\tConservG\tConservT\tQueryDel\tQueryIns\n" % (oligotype[oligoname][3], oligotype[oligoname][3]))
		        qPCResults[SetNbr][3][oligoname] =  {'Seq' : Seq, 'GC' : GC, 'LN' : LN, 'TM' : Tm, 'DG' : dG, 'ConsByPos' : {}, 'OligoScore' : 0, 'FreqArray' : []}
		        StartPos = int(qPCRData[qPCRset][oligotype[oligoname][4]])
		        EndPos = int(qPCRData[qPCRset][oligotype[oligoname][5]]) + 1
		        Score = 0
		        FreqArray = [None]*longest
		        if oligotype[oligoname][2] != "-":
			        for Posnbr in range(StartPos, EndPos):
				        qPCRsetOutFile.write("%sN%s\t%s\t%i\t" % (oligotype[oligoname][3],primerindex, Seq[primerindex -1 : primerindex], Posnbr))
				        Afreq, Cfreq, Gfreq, Tfreq, Delfreq, Insfreq, Nfreq = getFinalPosFreqs(SNPsDict, Posnbr, oligotype[oligoname][2])
				        qPCRsetOutFile.write("%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n" % (Afreq, Cfreq, Gfreq, Tfreq, Delfreq, Insfreq))
				        RefNT = RefGenomeSeq[Posnbr - 1]
				        qPCResults[SetNbr][3][oligoname]['ConsByPos'][primerindex] = [
				                Seq[primerindex -1 : primerindex], # Primer nucleotide at a given position
				                Posnbr, # Corresponding position in the reference genome
			                        Afreq, # Afreq at a given position
			                        Cfreq, # Cfreq at a given position
			                        Gfreq, # Gfreq at a given position
			                        Tfreq, # Tfreq at a given position
			                        Delfreq, # Delfreq at a given position
			                        Insfreq, # Insfreq at a given position
			                        Nfreq, # Degenerations frequency (N, K, M, D, etc.) at a given position
			                        RefNT  # Nt in Reference genome
				        ]
				        if primerindex > (len(Seq) - oligotype[oligoname][1]):
					        Consider = True			
				        if Consider == True:
					        weight = 1 / (EndPos - Posnbr)
					        Score += ConservDF.iloc[Posnbr-1]['Conservation'] * weight
				        FreqArray[primerindex -1] = sum(qPCResults[SetNbr][3][oligoname]['ConsByPos'][primerindex][col] for col in [value for value in NTColDict[qPCResults[SetNbr][3][oligoname]['ConsByPos'][primerindex][0]]])
				        primerindex += 1  
			        qPCResults[SetNbr][3][oligoname]['OligoScore'] = Score
			        qPCResults[SetNbr][3][oligoname]['FreqArray'] =  FreqArray   
		        else:
			        for Posnbr in reversed(range(StartPos, EndPos)):
				        qPCRsetOutFile.write("%sN%s\t%s\t%i\t" % (oligotype[oligoname][3],primerindex, Seq[primerindex -1 : primerindex], Posnbr))
				        Afreq, Cfreq, Gfreq, Tfreq, Delfreq, Insfreq, Nfreq = getFinalPosFreqs(SNPsDict, Posnbr, oligotype[oligoname][2])
				        qPCRsetOutFile.write("%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n" % (Afreq, Cfreq, Gfreq, Tfreq, Delfreq, Insfreq))
				        RefNT = rc_dna(RefGenomeSeq[Posnbr - 1])
				        qPCResults[SetNbr][3][oligoname]['ConsByPos'][primerindex] = [
				                Seq[primerindex -1 : primerindex], # Primer nucleotide at a given position
				                Posnbr, # Corresponding position in the reference genome
			                        Afreq, # Afreq at a given position
			                        Cfreq, # Cfreq at a given position
			                        Gfreq, # Gfreq at a given position
			                        Tfreq, # Tfreq at a given position
			                        Delfreq, # Delfreq at a given position
			                        Insfreq, # Insfreq at a given position
			                        Nfreq, # Degenerations frequency (N, K, M, D, etc.) at a given position
			                        RefNT # Nt in Reference genome
				        ] 
				        if primerindex > (len(Seq) - oligotype[oligoname][1]):
					        Consider = True			
				        if Consider == True:
					        weight = 1 / (EndPos - Posnbr)
					        Score += ConservDF.iloc[Posnbr-1]['Conservation'] * weight
				        FreqArray[primerindex -1] = sum(qPCResults[SetNbr][3][oligoname]['ConsByPos'][primerindex][col] for col in [value for value in NTColDict[qPCResults[SetNbr][3][oligoname]['ConsByPos'][primerindex][0]]])
				        primerindex += 1
			        qPCResults[SetNbr][3][oligoname]['OligoScore'] = Score
			        qPCResults[SetNbr][3][oligoname]['FreqArray'] =  FreqArray   
		                
	        TotalScore = 0
	        qPCRsetOutFile.write("#Final qPCR set conservation scores\n")
	        for oligoname in qPCResults[SetNbr][3]:
                        if TotalScore == 0:
                                qPCRsetOutFile.write("%sscore\t%0.6f" % (oligotype[oligoname][3], qPCResults[SetNbr][3][oligoname]['OligoScore']))
                                TotalScore += qPCResults[SetNbr][3][oligoname]['OligoScore']
                        else:
                                qPCRsetOutFile.write("\t%sscore\t%0.6f" % (oligotype[oligoname][3], qPCResults[SetNbr][3][oligoname]['OligoScore']))
                                TotalScore += qPCResults[SetNbr][3][oligoname]['OligoScore']
	        qPCResults[SetNbr][2] = TotalScore
	        qPCRsetOutFile.write("\nTotal score\t%0.6f" % (qPCResults[SetNbr][2]))

	        qPCRsetOutFile.close()
	        
	# Rank qPCR set by total score and report  
        qPCResults = collections.OrderedDict(sorted(qPCResults.items(), key=lambda t: t[1][2], reverse=True)) # only to order the dictionary
        qPCRRanking=open("qPCRRanking.txt", "w")
        print("\n\n*** REPORT ****\nqPCR set\tscore (Best scores first)")
        qPCRRanking.write("qPCR set\tscore\n")
        for key in qPCResults.keys():
	        qPCRRanking.write("%s\t%0.6f\n" % (qPCResults[key][0], qPCResults[key][2]))
	        print("%s\t%0.6f" % (qPCResults[key][0], qPCResults[key][2]))
        qPCRRanking.close()
        
        return qPCResults	       
        

### MAIN WORK ###
## Verbose section
# Starting screen follow-up
print("\nRunning ...")
LogFile=open("SNPs.log", "w") # for logging in future implementations

print("Getting annotations from GFF file %s ...                      " % (args.g), end="\r", flush=True)
CDSs = GFF2CDSs(GffFN) 
print("Getting annotations from GFF file %s ... Done!                " % (args.g))

# Import qPCRdata 
print("Reading qPCR information (file: %s) ...        " % (qPCRFN), end="\r", flush=True)
qPCRData = parse_qPCRdata(qPCRFN)
print("Reading qPCR information (file: %s) ... Done!       " % (qPCRFN))

### Creating some globally useful dictionaries and variables
# description of different types of oligos
oligotype = {"forward" : [9, getMinLen(qPCRData, 9), "+", "F", 4, 5, getMaxLen(qPCRData, 9)], # name : [sequence column in .tsv, min length, strand (+ ou -), single letter abbrev., start field, end field, max length]
     "reverse" : [15, getMinLen(qPCRData, 15), "-", "R", 11, 10, getMaxLen(qPCRData, 15)], 
     "probe" : [21, getMinLen(qPCRData, 21), "+", "P", 16, 17, getMaxLen(qPCRData, 21)] 
     }

#correspondance between primer bases and their frequencies array of qPCResults   
NTColDict = { 'A': [2], 'C': [3], 'G': [4], 'T': [5], 'R': [2, 4], 'Y': [3,5], 'S': [3,4], 'W': [2,5], 'K': [4,5], 'M': [2,3], 'B': [3,4,5], 'V': [2,3,4], 'D': [2,4,5], 'H': [2,3,5], 'N': [2,3,4,5] }

longest = max([oligotype[key][6] for key in oligotype.keys()]) # the length of the longest oligo

#Get reference genome data
Refrecords = getFastaData(RefGenoFN)
RefGenomeID = list(Refrecords[1].keys())[0]
RefGenomeSeq = Refrecords[1][RefGenomeID].replace("\n", "")
RefGenomeSize = len(RefGenomeSeq)
print("Assuming %s as reference (%i bases; Input file: %s)" % (RefGenomeID.rsplit(" ")[0], RefGenomeSize, args.i))

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

# Gather snp data using the filtered data frame
SNPsDict = getSNPs(GenomeDF, RefGenomeSeq)
del GenomeDF

# Retrieve a dataframe containing an array of conservation frequencies and an array of the corresponding positions in the refrecence genome
ConservDF = getConservArray(RefGenomeSeq, SNPsDict, FilteredNbrOfGeno, CDSs)  
# Retrieve and report global descriptive statistics and write figure
GlobalStatFile=open(("GlobalStatFile" + ".txt"), "w")
GlobalStatFile.write(str(ConservDF['Conservation'].describe(include='all')))
GlobalStatFile.write("\nmedian       %0.6f" % (ConservDF['Conservation'].median()))
GlobalStatFile.close()
print("\nDescriptive statistics of genomic data for %i genomes:" % (FilteredNbrOfGeno))
print (ConservDF['Conservation'].describe())
print("\nmedian       %0.6f\n" % (ConservDF['Conservation'].median()))

# Plot the conservation graph
fig = px.scatter(ConservDF, x="PosNbr", y="Conservation", color="Annotation(s)", labels={'PosNbr':'Position in the genome', 'Conservation':'Conservation frequency (0-1): ' }, hover_data=['PosNbr','Conservation', 'Annotation(s)', 'NTmutation(s)', 'AAmutation(s)' ], title="Conservation by position (Database date: %s) <br>Total number of genomes = %i" % (InfoDict['GDD'][1], FilteredNbrOfGeno))
fig.update_layout(
    hoverlabel_align = 'left')
fig.update_yaxes(range=[-0.1, 1.1])

# Write all conservation information to .tsv 
print("Writing conservation table: %s" % (OutFN), end="\r", flush=True)
with open(OutFN, 'w') as f:
	f.write("#START| **** INFORMATION SECTION ****\n")
	for code in InfoDict.keys():
		f.write("#%s|%s: %s\n" % (code, InfoDict[code][0], InfoDict[code][1]))
	f.write("#NFG| Number of filtered genomes: %s\n#AUF| Used filters: %s\n#DNA| DNA concentration: %0.2f\n#TMP| Temperature: %0.2f\n#NTP| dNTP concentration: %0.2f\n#MVI| Movalent ions concentration: %0.2f\n#DVI| Divalent ions concentration: %0.2f\n#END| **** INFORMATION SECTION ****\n" % (FilteredNbrOfGeno, ",".join(filters), DNAConc, Temp, dNTPConc, MVConc, DVConc))
	ConservDF.to_csv(f, index = False, sep="\t", mode='a')
print("Writing conservation table: %s - Done!" % (OutFN))

#fig.show()
fig.write_html("%s_ConservationByPos.html" % (InfoDict['GDD'][1]))

# Write output files, calculate qPCRset scores and return all colllected data
qPCResults = do_qPCR_analysis(qPCRData, ConservDF, RefGenomeSeq)

indexes = [qPCResults[SetNbr][0] for SetNbr in reversed(qPCResults.keys())]

# Create graph for visualization
for oligoname in oligotype.keys():
	FreqArray = []
	RefNT = []
	for SetNbr in reversed(qPCResults.keys()):
		FreqArray.append(qPCResults[SetNbr][3][oligoname]['FreqArray'])
		for primerindex in qPCResults[SetNbr][3][oligoname]['ConsByPos']:
			RefNT.append(qPCResults[SetNbr][3][oligoname]['ConsByPos'][primerindex][-1])

	### Draw heatmap creation
	fig = go.Figure(data=go.Heatmap(z=FreqArray, 
		x=[str(i-0.5) for i in range(1, len(FreqArray) + 1)],
		y=indexes,
		colorscale=[
		[0, "rgb(0, 0, 255)"],
		[0.8, "rgb(173,216,230)"],
		[0.99, "rgb(255,218,185)"],
		[1.0, "rgb(255, 0, 0)"]],
		colorbar={"title": 'Conservation frequency'},
		hovertemplate='qPCR set: %{y}<br>Primer position: %{x}<br>Primer NT frequency: %{z}',
		xgap = 1,
		ygap = 1,
		name = '',
		hoverongaps = False))
		
#	axis_template = dict(showgrid = False, 
##		gridwidth=1, 
##		gridcolor='White',
##		zeroline = True,
##		ticks = 'outside',
##		linecolor = 'grey', 
##		linewidth=1,
##		showticklabels = True,
##		mirror=False,
#		)
	xaxis=dict(
		title='Oligonucleotide base position',
		showgrid=False,
		tickmode = "array",
		showticklabels=True,
		ticks = 'outside',
		gridwidth=1, 
		gridcolor='White',
	)

	yaxis=dict(
		title='qPCR oligonucleotide set',
		showgrid=False,
		showticklabels=True,
	)

	fig.update_layout(title = "<b>Analysis of <i>{}</i> primers using <i>{}</i> dataset<br></b>".format(oligoname, InfoDict['GDD'][1]),
		  titlefont=dict(size =24, color='black', family='Arial, sans-serif'),
#		  margin = dict(t=100,r=100,b=100,l=100),
		  xaxis = xaxis,
		  yaxis = yaxis,
		  showlegend = True,
		  plot_bgcolor='rgb(233,233,233)',            
	#		  width = 800,
	#		  height = 1000,
		  autosize = True
	)

	#	fig.show()
	fig.write_html("%s_%s_primer_analysis.html" % (InfoDict['GDD'][1], oligoname))
	
					
LogFile.close()

print("\nJob finished !" )
quit()

