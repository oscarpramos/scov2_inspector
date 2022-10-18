#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### MODULES IMPORT ###
import argparse
import subprocess
import warnings
import gzip

warnings.filterwarnings("ignore")


### VARIABLES DEFINITION ###
nucmer = "nucmer"  # Path to nucmer executable - if in the system path then just the executable name is enough
show_snps = "show-snps" # Path to show-snps executable - if in the system path then just the executable name is enough
Nucmer_prefix = "Nucmer_out" 
SNPFileName = Nucmer_prefix + ".snp"
MetadataInfo = {}
Continent = OutFN = Country = Region = Conglomeration = Host = info = ""
DataDate = "YYYY-MM-DD"
TotalStartingGenomes = 0

# Help about how to use the script
parser = argparse.ArgumentParser(description='09a_getAllGenomesTable.py: Create a table containing all mutated genomes using multifasta file from GISAID and the matching metadata file. \n Usage: ./09a_getAllGenomesTable.py -i <MultiGenomeFile-GISAID(fasta)> -m <MetaData(tsv)> -r( <ReferenceFasta(fasta)> -g <GFFfile (gff)>\n', prefix_chars='-+')
parser.add_argument("-i", help="Genomes. The multifasta file download from GISAID that should be processed")
parser.add_argument("-m", help="Metadata provided by GISAID (tab separated)")
parser.add_argument("-r", help="The fasta file containing the reference genome")
parser.add_argument("-g", help="The gff file containing the annotation for the reference genome")
parser.add_argument("-o", help="Output file name for the compressed genome's table (default: GenomesTable.tsv.gz)")
parser.add_argument("-t", help="Number of threads (# of cores) to be used by nucmer during genome alignment phase (default = 2")

### FUNCTIONS ###
# Execute terminal commands and save the output to a logfile
def run(cmd, logfile):
        p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
        ret_code = p.wait()
        logfile.flush()
        return ret_code

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

# translate a DNA sequence to protein sequence
def translate_dna(sequence):
	proteinsequence = ''

	codontable = {
	    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	    }

	for n in range(0,len(sequence),3):
		if sequence[n:n+3] in codontable.keys():
		    proteinsequence += codontable[sequence[n:n+3]]
		else:
		    proteinsequence += "#" # non-expected codon (maybe codon composed of non standard DNA bases A, C, G, T or incomplete)
	sequence = ''
	return proteinsequence     

# get date format with fields separated by "-"	
def getDateFmt (string):
	DateElements = string.split("-")
	if len(DateElements) == 3 and DateElements[0].isdigit() == True and DateElements[1].isdigit() == True and DateElements[2].isdigit()  == True:
		Fmt = "%Y-%m-%d"
	elif len(DateElements) == 2 and DateElements[0].isdigit()  == True and DateElements[1].isdigit() == True:
		Fmt = "%Y-%m"
	elif DateElements[0].isdigit() == True:
		Fmt = "%Y"
	else:
		Fmt = "Invalid date format"

	return Fmt
	
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
        CDSHit = []
        for CDS in CDSs:
                if int(CDS[0]) <= int(now) <= int(CDS[1]):
                        if CDSHit == []:
                                CDSHit = [CDS]
                        else:
                                CDSHit.append(CDS)
        
        return CDSHit
	
# Get codon info from CDS with ribosome slippage (i.e, RDRP)
def get_RdRP_AAmut(Muta, CDS):
	Posnbr = int(Muta[1:-1])
	CDS_start = int(CDS[0])
	codonnbr = int(((Posnbr - CDS_start) / 3) + 1)
	SlippagePos = 13468
	SlippageCodonPos = 3
	SlippageLag = 1
	codon_info = []
	AAresult = ""
	if Posnbr < SlippagePos: # expected reading frame region
		codonnbr = int(((Posnbr - CDS_start) / 3) + 1)
		CodonPosEval = ((Posnbr + 1) - CDS_start ) / 3
		if CodonPosEval %1 == 0: 
			codon_info = [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr, 3] # codonnbr and positionnbr
		if 0.4 < CodonPosEval %1 < 0.8:
			codon_info = [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
		if 0 < CodonPosEval %1 < 0.4:
			codon_info = [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr, 1]
	elif  Posnbr == SlippagePos: # 2 codons are changed by mutations at this position
		codonnbr1 = int(((Posnbr - CDS_start) / 3) + 1)
		CodonPosEval1 = ((Posnbr + 1) - CDS_start ) / 3
		if CodonPosEval1 %1 == 0:
			codon_info = ['listoflists', [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr1, 3]] # codonnbr and positionnbr
		if 0.4 < CodonPosEval1 %1 < 0.8:
			codon_info = ['listoflists', [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr1, 2]]
		if 0 < CodonPosEval1 %1 < 0.4:
			codon_info = ['listoflists', [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr1, 1]]
		codonnbr2 = int(((Posnbr - (CDS_start - SlippageLag))/ 3) + 1)
		CodonPosEval2 = ((Posnbr + 1) - (CDS_start -1) ) / 3
		if CodonPosEval2 %1 == 0:
			codon_info.append([RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr2, 3]) # codonnbr and positionnbr
		if 0.4 < CodonPosEval2 %1 < 0.8:
			codon_info.append([RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr2, 2])
		if 0 < CodonPosEval2 %1 < 0.4:
			codon_info .append([RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr2, 1])	
	else: # shifted reding frame
		codonnbr = int(((Posnbr - (CDS_start - SlippageLag))/ 3) + 1)
		CodonPosEval = ((Posnbr + 1) - (CDS_start -1) ) / 3
		if CodonPosEval %1 == 0:
			codon_info = [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr, 3] # codonnbr and positionnbr
		if 0.4 < CodonPosEval %1 < 0.8:
			codon_info = [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
		if 0 < CodonPosEval %1 < 0.4:
			codon_info = [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr, 1]
		
	# generate new codon from mutation, translate and report the AA mutation.
	if codon_info[0] != 'listoflists':
		new_codon = list(codon_info[0])
		new_codon[codon_info[2] - 1] = Muta[-1]
		new_codon = "".join(new_codon)
		AAresult = "%s%i%s" % (translate_dna(codon_info[0]), codon_info[1], translate_dna(new_codon))
	else:
		for muta in codon_info[1:]:
			new_codon = list(muta[0])
			new_codon[muta[2] - 1] = Muta[-1]
			new_codon = "".join(new_codon)
			AAmut = "%s%i%s" % (translate_dna(muta[0]), muta[1], translate_dna(new_codon))
			if AAresult == "":
				AAresult = AAmut
			else:
				AAresult = AAresult + '&' + AAmut	

	return AAresult

# Translate the consequence of NT mutations to aa (general CDSs)
def get_AAmut(Muta, CDS):
        Posnbr = int(Muta[1:-1])
        CDS_start = int(CDS[0])
        codon_info = []
        codonnbr = int(((Posnbr - CDS_start) / 3) + 1)
        CodonPosEval = ((Posnbr + 1) - CDS_start ) / 3
        if CodonPosEval %1 == 0:
                codon_info = [RefGenomeSeq[Posnbr -3 : Posnbr], codonnbr, 3] # codon, codonnbr and positionnbr
        if 0.4 < CodonPosEval %1 < 0.8:
                codon_info = [RefGenomeSeq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
        if 0 < CodonPosEval %1 < 0.4:
                codon_info = [RefGenomeSeq[Posnbr -1 : Posnbr + 2], codonnbr, 1]
        
        # generate new codon from mutation, translate and report the AA mutation. 
        new_codon = list(codon_info[0])
        new_codon[codon_info[2] - 1] = Muta[-1]
        new_codon = "".join(new_codon)
        AAresult = "%s%i%s" % (translate_dna(codon_info[0]), codon_info[1], translate_dna(new_codon))
        
        return AAresult
        
# return full mutation information with annotation	
def get_ProtMutas(NtMuta, CDSs):
	ProtMuta = ""
	RefPos = NtMuta[1:-1]
	AAResult = ""
	Annot = ""
	AnnotCDSs = IFCDS(CDSs, int(RefPos)) # get annotation associated to the mutation

	if AnnotCDSs != []:
		for CDS in AnnotCDSs:	
			if CDS[2] == "nsp12(RdRP) ":
				AAResult = get_RdRP_AAmut(NtMuta, CDS)
				Annot = CDS[2].strip()
			elif CDS[2] == "5'UTR " or CDS[2] == "3'UTR ":
				Annot = CDS[2].strip()
				AAResult = "-"
			else: 
				AAResult = get_AAmut(NtMuta, CDS)  # get AAmutation associated to each annotation
				Annot = CDS[2].strip()
	else:
		Annot = "UTR"
		AAResult = "-"

#	print(AAResult)
	ProtMuta = AAResult + ";" + Annot

	return ProtMuta	

# Load FastaDict from fasta file	
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
			print("  Reading genomes from file %s ...                    " % (FastaFileName), end="\r", flush=True)
			PercLast = PercNow
	print("  Reading genomes from file %s ... Done!                    " % (FastaFileName))
		
	return a, FastaDict

#Count the number of degenerated bases for each sequence in a multifasta file
def getDegPos(FastaFileName):
        DegPDict = {}
        FastaRecords = getFastaData(FastaFileName)
        for header in FastaRecords[1].keys():
                counter = 0

                for char in FastaRecords[1][header]:
                        if char not in ['A','C','G','T']:
                                counter += 1
                DegPDict[header] = counter
                
        del FastaRecords
        
        return DegPDict
                              
# Retrieve metadata information
def getMetadataInfo(MetaFN):
	MetaInfo = {}
	extension = MetaFN.split(".")[-1]
	if  extension == "gz" or extension == "gzip":
		Metadata = get_gzfile_content(MetaFN)
	else:
		Metadata = get_file_content(MetaFN)
		
	for line in Metadata[1:]:
		LineFields = line.split("\t")
		MetaVirus = LineFields[0].strip()
		MetaAccID = LineFields[2].strip()
		MetaCollDate = LineFields[4].strip()
		MetaConti = LineFields[9].strip()
		MetaCountry = LineFields[10].strip()
		MetaRegion = LineFields[11].strip()
		MetaSeqLen = LineFields[13].strip()
		MetaHost = LineFields[14].strip()
		MetaPangLine = LineFields[18].strip()
		if MetaPangLine == "" or MetaPangLine == None or MetaPangLine == " ":
		        MetaPangLine = "Not known"
		MetaSubDate = LineFields[26].strip()
		if MetaVirus not in MetaInfo.keys():
			MetaInfo[MetaVirus] = { 
			'MetaAccID' : MetaAccID,
			'MetaCollDate' : MetaCollDate,
			'MetaSubDate' : MetaSubDate,
			'MetaConti' : MetaConti,
			'MetaCountry' : MetaCountry,
			'MetaRegion' : MetaRegion,
			'MetaHost' : MetaHost,
			'MetaPangLine' : MetaPangLine,
			'MetaSeqLen' : MetaSeqLen,
			'MetaHost' : MetaHost,
			'MetaPangLine' : MetaPangLine,
			'MetaSubDate' : MetaSubDate
		}
		else:
			print("Duplicated virus : %s " % (MetaVirus))

	return MetaInfo
	
# Get mutations reported in show-aligns output
def get_SNPs(fileName, CDSs):
	SNPdata_00 = {}
	try:
		SNPinfo = get_file_content(fileName)[4:] # from the 5th line to the end
		a = 0
		PercLast = 0
		TotalRecs = len(SNPinfo)
		for line in SNPinfo:
			Fields = line.rsplit("\t")
			QueryNuc = Fields[2]
			if QueryNuc in ['A', 'C', 'G', 'T', '.']: # Exclude non-determined nucleotides (N for instance) 
				RefNuc = Fields[1]
				RefPos = Fields[0]
				
				QueryGenomeID = Fields[11].rsplit("|")[0]
				NTmuta = RefNuc + RefPos + QueryNuc
				AAmuta = get_ProtMutas(NTmuta, CDSs)
				FullMutaDesc = NTmuta + ";" + AAmuta 
				if QueryGenomeID in SNPdata_00.keys():	
					SNPdata_00[QueryGenomeID].append(FullMutaDesc)
				else:
					SNPdata_00[QueryGenomeID] = [FullMutaDesc]
				a += 1
				PercNow = int((a * 100)/TotalRecs)
				if PercLast != PercNow:
					print("  Gathering SNP data ... %i %s                                                   " % (PercNow, "%"), end="\r", flush=True)
					PercLast = PercNow
			else:
				continue
		print("  Gathering SNP data ... Done!      ") 
		return SNPdata_00
		
	except ValueError:
                return "Could not get SNPs!"
	
def doMummer(MultiFastaFile, CDSs):
	## Mummer section ##
	# run Nucmer using the corrected fasta file 
	cmd = ("%s -t %i -p %s %s %s  2>/dev/null" % (nucmer, Nthreads, Nucmer_prefix, args.r, MultiFastaFile)) # Redirect errors to /dev/null; Nucmer version 4rc1
	print("  Running Nucmer ...             ", end="\r", flush=True)
	run(cmd, LogFile)
	print("  Running Nucmer ... Done !       ")
	# run show-snps and redirect the output to SNPfile
	cmd = ("%s %s %s 2>/dev/null" % (show_snps, "-q -T", Nucmer_prefix + ".delta")) # Redirect errors to /dev/null
	SNPFile=open(SNPFileName, "w")
	print("  Running show-snps ...           ", end="\r", flush=True)
	run(cmd, SNPFile)
	SNPFile.close()
	print("  Running show-snps ... Done!      ")
	# read SNPfile and put information into a dictionary and do some work
	SNPdata = get_SNPs(SNPFileName, CDSs)

	return SNPdata

def writeAllSNPdata(AllVarData, MetaDataInfo, DegPDict, RefGenomeID, RefGenomeSize, TotalStartingGenomes, DataDate, gzipped):
	if gzipped == True:
		AllSNPGenFile = gzip.open("IndividualGenomeList.tsv.gz", 'wt')
	else:
		AllSNPGenFile=open("IndividualGenomeList.tsv", "w")
	AllSNPGenFile.write("#START| **** INFORMATION SECTION ****\n#GDD| GISAID database date: %s\n#DFN| Input file: %s\n#NoG| Number of genomes: %i\n#RID| Reference Genome: %s\n#RLE| Number of bases: %i\n#RFN| Reference genome file: %s\n#AFN| Reference annotation file: %s\n#MFN| Metadata file: %s\n#END| **** INFORMATION SECTION ****\n" % (DataDate, args.i, TotalStartingGenomes, RefGenomeID.rsplit(" ")[0], RefGenomeSize, args.r, args.g, args.m))
	AllSNPGenFile.write("Name\tGISAID_ID\t#Mutations\tPointMutations\tCollection_Date\tSubmission_Date\tContinent\tCountry\tRegion\tHost\tClade\tSequence_length\t#Ns")
	for Virus in AllVarData.keys():
		AllSNPGenFile.write("\n%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i" % (Virus, MetaDataInfo[Virus]['MetaAccID'], len(AllVarData[Virus]), "|".join(AllVarData[Virus]), 
		MetaDataInfo[Virus]['MetaCollDate'], MetaDataInfo[Virus]['MetaSubDate'], MetaDataInfo[Virus]['MetaConti'], MetaDataInfo[Virus]['MetaCountry'], MetaDataInfo[Virus]['MetaRegion'],
		MetaDataInfo[Virus]['MetaHost'], MetaDataInfo[Virus]['MetaPangLine'], MetaDataInfo[Virus]['MetaSeqLen'], DegPDict[Virus]))


### INITIALIZING TASKS ###
print("")
#actually parse command line options
args = parser.parse_args()
if args.i is not None:
	MultiFastaFN = args.i
	SeqDate = MultiFastaFN.rsplit("/")[-1].rsplit("_")[1]
	SeqDateFmt = getDateFmt(SeqDate)
	if SeqDateFmt != "Invalid date format":
		 DataDate = SeqDate
		 del SeqDate, SeqDateFmt
else:
	print("  Please provide genomes in fasta format provided by GISAID (-i)")
	quit()

if args.m is not None:
	MetaDataFN = args.m
else:
	print("  Please provide metadata provided by GISAID (tab separated)  (-m)")
	quit()
	
if args.r is not None:
	RefGenome = args.r
else:
	print("  Please provide a reference genome (-r)")
	quit()
	
if args.g is not None:
	GffFN = args.g
else:
	print("  Please provide gff annotation file (-g)")	
	quit()
	
if args.o is not None:
	OutFN = args.o
	extension = args.o.rsplit(".")[-1]
	if extension != "gz" and  extension != "gzip":
	      OutFN = OutFN + '.gz'  
else:
        OutFN = "IndividualGenomesTable.tsv.gz"
        
if args.t is not None:
	Nthreads = int(args.t)
else:
        Nthreads = 2

#Get number of "N"s for each sequence
print("  Count the number of degenerated bases in file %s ...                      " % (args.i), end="\r", flush=True)
DegPDict = getDegPos(args.i)
print("  Count the number of degenerated bases in file %s ... Done!                      " % (args.i))

#Get reference genome data
Refrecords = getFastaData(args.r)
RefGenomeID = list(Refrecords[1].keys())[0]
RefGenomeSeq = Refrecords[1][RefGenomeID].replace("\n", "")
RefGenomeSize = len(RefGenomeSeq)
print("  Assuming %s as reference (%i bases; Input file: %s)" % (RefGenomeID.rsplit(" ")[0], RefGenomeSize, args.i))

# Evalutate the number of starting genomes
TotalStartingGenomes = int(subprocess.check_output(["grep", "-c", ">", args.i]).decode("utf-8")) 
print("  %s starting genomes found!" % TotalStartingGenomes)

print("  Getting annotations from GFF file %s ...                      " % (args.g), end="\r", flush=True)
CDSs = GFF2CDSs(GffFN) 
print("  Getting annotations from GFF file %s ... Done!                " % (args.g)) 

	
### MAIN WORK ###
## Verbose section
# Starting screen follow-up
print("\nRunning ...")
Phase = "Global tasks"
print("%s" % (Phase))
LogFile=open("IndivGen.log", "w")

# Run mummer and get the list of SNPs / genome
AllSNPdata = doMummer(args.i, CDSs)

# Gather metadata information
print("  Getting metadata from file %s  ...                    " % (args.m), end="\r", flush=True)
MetaInfo = getMetadataInfo(args.m)
print("  Getting metadata from file %s  ... Done!              " % (args.m))

# Write genome information
gzipped = True
writeAllSNPdata(AllSNPdata, MetaInfo, DegPDict, RefGenomeID, RefGenomeSize, TotalStartingGenomes, DataDate, gzipped)

LogFile.close()

print("Job finished !" )
quit()

