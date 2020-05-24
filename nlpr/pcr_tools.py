#pcr_tools.py
# THIS PROGRAM WAS ORIGINALLY WRITTEN BY kmb* sept 20, 2012
# BUT THE FILEPATHS WERE UPDATED TO RUN ON MY MACHINE (JGrembi) ON  29 July 2014
primer3_path = '/usr/bin/' #/usr/bin/primer3_core

def hello_pcr_tools():
	print "## \nYou are using a version pcr_tools.py last updated on Sept 20, 2012 \n##"

def file_to_list(sys_index = 1,line_index = 0):
	L = list()
	import sys
	fn = sys.argv[sys_index]
	fh = open(fn , 'r')
	for line in fh:
		L.append(line.strip().split()[line_index])
	return L

def generate_seq_dictionary(x):
	# arg: a <filename.fna>
	# returns: a dictionary <record_dict>
	# USE ME IF YOU WANT TO GENERATE A DICTIONARY OF SEQUENCE
	# example: master_dictionary = generate_seq_dictionary(sys.argv[1])
	from Bio import SeqIO
	fh = open(x, "rU")
	record_dict = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
	fh.close()
	print "## \nYou generated a sequence dictionary with %i entries \n##"%(len(record_dict.keys()))
	return record_dict


def write_subfasta(L,D,filename):
	# args: 
		# (1) a List of the accession you want to create
		# (2) Dictionary to search
		# (3) filename to write output file  
	# returns 
		# write the output file 
	# USE ME WITH THE ABOVE generate_seq_dictionary TO MAKE A SUB FASTA FILE FROM THE MASTER DICTIONARY
		# example: 
			# master_dictionary = generate_seq_dictionary(sys.argv[1])
			# write_subfasta(list_of_accessions, master_dictionary, 'output.fna')
	fh = open(filename, "w") 
	for l in L:
		#print "> %s \n%s\n" %(D[l].id,str(D[l].seq))
		try:
			fh.write("> %s \n%s\n" %(D[l].id,str(D[l].seq)))
		except KeyError:
			continue
	fh.close()
	print "## \nYou wrote a subfasta with with %i sequence(s)\n##"%(len(L))

def write_full_subfasta(L,D,filename):
	# args: 
		# (1) a List of the accession you want to create
		# (2) Dictionary to search
		# (3) filename to write output file  
	# returns 
		# write the output file 
	# USE ME WITH THE ABOVE generate_seq_dictionary TO MAKE A SUB FASTA FILE FROM THE MASTER DICTIONARY
		# example: 
			# master_dictionary = generate_seq_dictionary(sys.argv[1])
			# write_subfasta(list_of_accessions, master_dictionary, 'output.fna')
	fh = open(filename, "w") 
	for l in L:
		#print "> %s \n%s\n" %(D[l].id,str(D[l].seq))
		try:
			fh.write(">%s \n%s\n" %(D[l].description,str(D[l].seq)))
		except KeyError:
			continue
	fh.close()
	print "## \nYou wrote a subfasta with with %i sequence(s)\n##"%(len(L))


def fna_to_pr3_input(In, Out):
	# For now the parameters are hardcoded, but this takes a fasta file and outputs .pr3in file
	import sys
	from Bio import SeqIO
	file_handle = open(In, 'r')
	output_handle = open(Out, "w")
	#This block writes the parameters to be use by primer3. Version 1.4
	output_handle.write("PRIMER_NUM_RETURN=500\n") 
	output_handle.write("PRIMER_MIN_TM=54\n")
	output_handle.write("PRIMER_OPT_TM=55\n")
	output_handle.write("PRIMER_MAX_TM=56\n")
	#output_handle.write("PRIMER_OPT_SIZE=17\n")
	output_handle.write("PRIMER_MIN_SIZE=15\n")
	output_handle.write("PRIMER_MAX_SIZE=26\n")
	output_handle.write("PRIMER_NUM_NS_ACCEPTED=1\n")
	output_handle.write("PRIMER_PRODUCT_SIZE_RANGE=75-200\n")
	output_handle.write("PRIMER_GC_CLAMP=0\n")
	output_handle.write("PRIMER_FILE_FLAG=1\n")
	output_handle.write("PRIMER_EXPLAIN_FLAG=1\n")
	#This block writes the sequences in sys.argv[1] to the primer3input file
	counter = 0 
	for seq_record in SeqIO.parse(file_handle, 'fasta'):
		output_handle.write("PRIMER_SEQUENCE_ID=" + str(seq_record.id) +"\n") 
		output_handle.write("SEQUENCE=" + str(seq_record.seq) + "\n")
		output_handle.write("="+ "\n")
		counter = counter + 1
	print "## \nYou generated a primer3 input file with %i sequence(s)\n##"%(counter)


def fna_to_pr3_v2_3_4_input(In, Out, number_of_primers_to_generate):
	# For now the parameters are hardcoded, but this takes a fasta file and outputs .pr3in file
	import sys
	from Bio import SeqIO
	file_handle = open(In, 'r')
	output_handle = open(Out, "w")
	#This block writes the parameters to be use by primer3. Version 1.4
	output_handle.write("PRIMER_NUM_RETURN=%i\n"%(number_of_primers_to_generate)) 
	output_handle.write("PRIMER_MIN_TM=59\n")
	output_handle.write("PRIMER_OPT_TM=60\n")
	output_handle.write("PRIMER_MAX_TM=61\n")
	output_handle.write("PRIMER_MIN_SIZE=15\n")
	output_handle.write("PRIMER_MAX_SIZE=28\n")
	output_handle.write("PRIMER_NUM_NS_ACCEPTED=1\n")
	output_handle.write("PRIMER_PRODUCT_SIZE_RANGE=45-200\n") # CHANGED FOR tRNA study Feb 10, 2014
	output_handle.write("PRIMER_GC_CLAMP=0\n")
	output_handle.write("PRIMER_FILE_FLAG=1\n")
	output_handle.write("PRIMER_EXPLAIN_FLAG=1\n")
	output_handle.write('PRIMER_TM_FORMULA=1\n')
	output_handle.write('PRIMER_SALT_CORRECTIONS=1\n')
	output_handle.write('PRIMER_THERMODYNAMIC_ALIGNMENT=1\n')
	output_handle.write('PRIMER_SALT_DIVALENT=3\n')
	output_handle.write('PRIMER_DNTP_CONC=0.6\n')
	output_handle.write('PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0\n')
    output_handle.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=' + primer3_path + 'primer3_config/\n')
    # output_handle.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/koshlan/primer3-2.3.4/src/primer3_config/\n') #!UPDATE ME!#
	#output_handle.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/JGrembi/primer3-2.3.4/src/primer3_config/\n')
	#This block writes the sequences in sys.argv[1] to the primer3input file
	
	counter = 0 
	for seq_record in SeqIO.parse(file_handle, 'fasta'):
		output_handle.write("SEQUENCE_ID=" + str(seq_record.id) +"\n") 
		output_handle.write("SEQUENCE_TEMPLATE=" + str(seq_record.seq) + "\n")
		output_handle.write("="+ "\n")
		counter = counter + 1
	print "## \nYou generated a primer3v2.3.4 input file with %i sequence(s)\n##"%(counter)


def write_check_p3(p1,p2,Out):
	oh = open(Out, 'w')
	oh.write("PRIMER_TASK=check_primers\n")
	oh.write("SEQUENCE_PRIMER=%s\n"%(p1))
	oh.write("SEQUENCE_PRIMER_REVCOMP=%s\n"%(p2))
	oh.write("PRIMER_MIN_TM=55\n")
	oh.write("PRIMER_MAX_TM=70\n")
	oh.write("PRIMER_FILE_FLAG=1\n")
	oh.write("PRIMER_EXPLAIN_FLAG=1\n")
	oh.write("PRIMER_TM_FORMULA=1\n")
	oh.write("PRIMER_SALT_CORRECTIONS=1\n")
	oh.write("PRIMER_THERMODYNAMIC_ALIGNMENT=1\n")
	oh.write("PRIMER_SALT_DIVALENT=3\n")
	oh.write("PRIMER_DNTP_CONC=0.6\n")
	oh.write("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0\n")
	oh.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + primer3_path + "primer3_config/\n")
    # oh.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/JGrembi/primer3-2.3.4/src/primer3_config/\n")#!UPDATE ME!#
    # oh.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/koshlan/primer3-2.3.4/src/primer3_config/\n")
	oh.write("=\n")

def pr3_modern(In, Out, pr3_path):
	import os
	os.system(pr3_path + "primer3_core <" + In + ">" + Out)
    
def pr3(In, Out):
	import os
	os.system(primer3_path + "primer3_core <" + In + ">" + Out)
    #os.system("/Users/JGrembi/primer3-2.3.4/src/primer3_core <" + In + ">" + Out)
    # os.system("/Users/koshlan/primer3-2.3.4/src/primer3_core <" + In + ">" + Out) #!UPDATE ME!#
    # JULY 24, 2014 (I MADE THIS CHAGEN SO JESS COULD RUN ON HER MACHINE)#os.system("/Users/koshlan/primer3-2.3.4/src/primer3_core <" + In + ">" + Out)
	#os.system("primer3_core <" + In + ">" + Out)


def parse_pr3out(In, Out):
	# This was used from January 15, call_primer3_v3.py should work fine
	primer3output = open(In, 'r')
	output_handle = open(Out, 'w')
	primers_dictionary = {}
	PRIMER_PAIR = '0'
		#######################################################
		#### II.i if statements to give each sequence an ID ###
		#######################################################
	for line in primer3output:
		# if the line starts with the sequence identifier (typically first line), we store that information as a variable
		if line.startswith("PRIMER_SEQUENCE_ID"):
			#define the parent sequence name as the second element after line is split it into 2 elements on either side of the "=" symbol and strip away the newline
			PRIMER_SEQUENCE_ID = line.split("=")[1].strip("\n") 
			continue 
		#as we parse each pair begins with a PRIMER_PAIR_PENALTY_ Designation
		# for the first case 
		if line.startswith("PRIMER_PAIR_PENALTY") and line.split('=')[0] == "PRIMER_PAIR_PENALTY":
			PRIMER_PAIR = '0'
			PRIMER_PAIR_PENALTY = line.split('=')[1].strip('\n')
			continue					
		if line.startswith("PRIMER_PAIR_PENALTY_"): 
			# strip away the Primer_Pair_Penalty portion of the string leaving just the number which will be use to identify the pair
			try:
				PRIMER_PAIR = line.split("=")[0]
				PRIMER_PAIR = PRIMER_PAIR.split("PRIMER_PAIR_PENALTY_")[1]
				PRIMER_PAIR_PENALTY = line.split('=')[1].strip('\n')
			except IndexError:
				continue
		#######################################################
		#### II.ii if statements to get sequence and make dictionary
		#######################################################
		if "PRIMER_LEFT" in line and "SEQUENCE" in line:
			PRIMER_LEFT_SEQUENCE = line.split("=")[1].strip('\n')
			primers_dictionary[PRIMER_SEQUENCE_ID + "_" + PRIMER_PAIR] = [PRIMER_LEFT_SEQUENCE] 
		if "PRIMER_RIGHT" in line and "SEQUENCE" in line:
			PRIMER_RIGHT_SEQUENCE = line.split("=")[1].strip('\n')
			primers_dictionary[PRIMER_SEQUENCE_ID + "_" + PRIMER_PAIR].append(PRIMER_RIGHT_SEQUENCE) 
		#######################################################
		#### II.iii if statements to parse Tm, primer_start_positions, and penalty_score and PRINT
		#######################################################
		#SPECIAL if Statments for the first case of primer position
		if line.split("=")[0] == "PRIMER_LEFT":
			LEFT_START = line.split("=")[1].split(',')[0]	
		if line.split("=")[0] == "PRIMER_RIGHT":
			RIGHT_START = line.split("=")[1].split(',')[0]
		#SPECIAL if Statments for 2nd through nth case for primer position
		if line.split('=')[0] == "PRIMER_LEFT_" + PRIMER_PAIR:
			LEFT_START = line.split("=")[1].split(',')[0]	
		if line.split('=')[0] == "PRIMER_RIGHT_" + PRIMER_PAIR:
			RIGHT_START = line.split("=")[1].split(',')[0]		
		#if Statments for 1st through nth case for TM
		if "PRIMER_LEFT" in line and 'TM' in line:
				LEFT_TM = line.split("=")[1].strip('\n')
		if "PRIMER_RIGHT" in line and 'TM' in line:
				RIGHT_TM = line.split("=")[1].strip('\n')
		if "PRIMER_PRODUCT_SIZE" in line and not("RANGE" in line):
			PRIMER_PRODUCT_SIZE = line.split("=")[1].strip('\n')
		#######################################################
		#### II.iii PRINT and WRITE OUTPUT
		#######################################################
			#print ">" + PRIMER_SEQUENCE_ID + "_" + PRIMER_PAIR + "\t"+ PRIMER_PAIR_PENALTY + "\t" + PRIMER_LEFT_SEQUENCE + "\t" + PRIMER_RIGHT_SEQUENCE +"\t"+ LEFT_TM + "\t" + RIGHT_TM + "\t" + LEFT_START + "\t" + RIGHT_START + "\t"+ PRIMER_PRODUCT_SIZE   
			output_handle.write(">" + PRIMER_SEQUENCE_ID + "_" + PRIMER_PAIR + "\t"+ PRIMER_PAIR_PENALTY + "\t" + PRIMER_LEFT_SEQUENCE + "\t" + PRIMER_RIGHT_SEQUENCE +"\t"+ LEFT_TM + "\t" + RIGHT_TM + "\t" + LEFT_START + "\t" + RIGHT_START + "\t"+ PRIMER_PRODUCT_SIZE + "\n")


def parse_pr3_v2_3_4_out(In, Out):
	# This was used from January 15, call_primer3_v3.py should work fine
	primer3output = open(In, 'r')
	output_handle = open(Out, 'w')
	primers_dictionary = {}
	PRIMER_PAIR = '0'
		#######################################################
		#### II.i if statements to give each sequence an ID ###
		#######################################################
	import re
	for line in primer3output:
		# if the line starts with the sequence identifier (typically first line), we store that information as a variable
		if line.startswith("SEQUENCE_ID"):
			#define the parent sequence name as the second element after line is split it into 2 elements on either side of the "=" symbol and strip away the newline
			PRIMER_SEQUENCE_ID = line.split("=")[1].strip("\n") 
			continue 
		if re.search('PRIMER_PAIR_([\d]+)_PENALTY', line): # SEARCH FOR "PRIMER_PAIR_1_PENALTY"
			matchObj = re.search('PRIMER_PAIR_([\d]+)_PENALTY', line)
			try:
				PRIMER_PAIR = matchObj.group(1)
				PRIMER_PAIR_PENALTY = line.split('=')[1].strip('\n')
			except IndexError:
				continue
		#######################################################
		#### II.ii if statements to get sequence and make dictionary
		#######################################################
		if re.search('PRIMER_LEFT_([\d]+)_SEQUENCE', line):
			PRIMER_LEFT_SEQUENCE = line.split("=")[1].strip('\n')
			continue
		if re.search('PRIMER_RIGHT_([\d]+)_SEQUENCE', line):
			PRIMER_RIGHT_SEQUENCE = line.split("=")[1].strip('\n')
			continue
		if re.search('PRIMER_LEFT_([\d]+)=', line):
			LEFT_START = line.split("=")[1].split(',')[0]
			continue
		if re.search('PRIMER_RIGHT_([\d]+)=', line):
			RIGHT_START = line.split("=")[1].split(',')[0]
			continue
		if re.search('PRIMER_LEFT_([\d]+)_TM', line):
			LEFT_TM = line.split("=")[1].strip('\n')
			continue
		if re.search('PRIMER_RIGHT_([\d]+)_TM', line):
			RIGHT_TM = line.split("=")[1].strip('\n')
			continue
		if re.search('PRIMER_PAIR_([\d]+)_PRODUCT_SIZE', line):
			PRIMER_PRODUCT_SIZE = line.split("=")[1].strip('\n')
			output_handle.write(">" + PRIMER_SEQUENCE_ID + "_" + PRIMER_PAIR + "\t"+ PRIMER_PAIR_PENALTY + "\t" + PRIMER_LEFT_SEQUENCE + "\t" + PRIMER_RIGHT_SEQUENCE +"\t"+ LEFT_TM + "\t" + RIGHT_TM + "\t" + LEFT_START + "\t" + RIGHT_START + "\t"+ PRIMER_PRODUCT_SIZE + "\n")
			continue
			
def write_fuznuc_pattern(p1,p2, Out, mismatches = 4):
	#args: 
		# p1 (string) first primer sequence
		# p2 (string) second primer sequence 
		# mismatches (int) number of acceptable mismatches per primer (DEFAULT = 4)
		# Out (string) name of temporary file to write the pattern file
	# results: writes a temporary file 
		# 	>pat1 <mismatch=4> 
		# 	GTTAGTCCCTTGGTGGCT       
		# 	>pat2 <mismatch=4> 
		# 	CGGGTCTAAAGCCTTTTC 
	fh = open(Out, 'w')
	fh.write(">pat1 <mismatch=%i>\n%s\n>pat2 <mismatch=%i>\n%s\n"%(mismatches,p1,mismatches,p2))



def execute_fuznuc(pattern_file,target_seqs, Out):
	# Args:
		# pattern_file (string): name of the temporary pattern file .pat
		# target (string): name of the sequences file .fna
		# Out (sting): name of the temporary output file
	#NOTE: ENSURE THE FOLLOWING IS INSTALLED in ~/EMBOSS-6.4.0/emboss/fuzznuc
	#NOTE: http://emboss.sourceforge.net/docs/themes/ReportFormats.html
	import os
	os.system('~/EMBOSS-6.5.7/emboss/fuzznuc %s %s -pattern @%s -complement Y -rformat excel'%(target_seqs, Out, pattern_file))
	#print "## \nFrom %s and %s You generated a fuznuc output file: %s\n##"%(target_seqs, pattern_file,Out)
	#EXAMPLE:
	#os.system('/Users/koshlan/EMBOSS-6.4.0/emboss/fuzznuc fuzzIn.fasta fuzzOut.pat -pattern @oligo_nucseq.pat -complement Y')


def parse_fzout(In):
	# args: In (string) f
	# Returns Dictionary D[SeqName] {   {   'AB301952.1': {   'negative': '1', 'positive': '0', 'total': 0},	
	fh = open(In, 'r')
	D={}
	for line in fh:
		# SKIP HEADER LINES
		if line.startswith("SeqName"):
			continue
		else:
			line = line.strip()
			L = line.split()
			(SeqName,Start,End,Score,Strand,Pattern_Seq,Mismatch) = line.split()
			Pattern = Pattern_Seq.split(':')[0] # KMB JAN 23, 2013, fixed the fact that pat:seq 
			Seq = Pattern_Seq.split(':')[1] 
			Mismatch = Mismatch.replace('.', '0') # Changed out the . for a zeor mismatch
			if SeqName not in D.keys():
				D[SeqName] = dict()
				D[SeqName]['pat1'] = 99
				D[SeqName]['pat2'] = 99
			if Pattern == 'pat1': # KMB JAN 27, 2013 CHANGED to 'pat1:' to 'pat1'
				D[SeqName]['pat1'] = int(Mismatch)
			elif Pattern == "pat2": # KMB JAN 27, 2013 CHANGED to 'pat1:' to 'pat1'
				D[SeqName]['pat2'] = int(Mismatch)
	L = D.keys()
	for l in L:
		try: 
			X = int(D[l]['pat1']) + int(D[l]['pat2'])
		except KeyError:
			D[l]['total'] = 99
		else: 
			D[l]['total'] = X
	#print "## \nThe function <parse_fzout> recieved %s, and returned a dictionary"%(In)
	# import pprint
	# pp = pprint.PrettyPrinter(indent=4)
	# pp.pprint(D)
	return D


def hit_percentage_calc(L,D,acceptable_mismatches):
	hit_count = 0
	for l in L:
		try:
			(l in D.keys()) & (D[l]['total'] < acceptable_mismatches)
		except KeyError:
			continue
		else:
		 	if (l in D.keys()) & (D[l]['total'] < acceptable_mismatches):
		 		hit_count = hit_count + 1
	
	hit_percentage = round(float(hit_count)/float(len(L)), 3) # (1)
	return hit_percentage


def hit_percentage_calc_plus_discovery(L,D,acceptable_mismatches):
	hit_count = 0
	Lout = []
	for l in L:
		try:
			(l in D.keys()) & (D[l]['total'] < acceptable_mismatches)
		except KeyError:
			continue
		else:
		 	if (l in D.keys()) & (D[l]['total'] < acceptable_mismatches):
		 		Lout.append(l)
				hit_count = hit_count + 1
	hit_percentage = round(float(hit_count)/float(len(L)),3) # (1)
	#print hit_percentage
	#print Loutl
	#print L
	#print D
	return {'hit_percentage' : hit_percentage, 'hits': Lout}
	

def hit_percentage_calc_plus_number_of_mismatch(L,D,acceptable_mismatches):
	hit_count = 0
	Lout = []
	Lout_mismatch = []
	for l in L:
		try:
			(l in D.keys()) & (D[l]['total'] < acceptable_mismatches)
		except KeyError:
			continue
		else:
		 	if (l in D.keys()) & (D[l]['total'] < acceptable_mismatches):
		 		Lout.append(l)
				Lout_mismatch.append(D[l]['total'])
				hit_count = hit_count + 1
	hit_percentage = round(float(hit_count)/float(len(L)),3) # (1)
	#print hit_percentage
	#print Loutl
	#print L
	#print D
	return {'hit_percentage' : hit_percentage, 'hits': Lout, 'number_mismatched': Lout_mismatch }


def sort_single_capture_primers(In,Out,number_of_primers_to_review):
	# this finds the primers that best capture the whole cluster, return the 100 best.
	import os
	os.system("sort -k 10,10nr -k 2,2n %s > tempsort.txt" %(In)) # SORT COMMAND, SORT FIRST ON THE 10th column numerically, next sort on the 2nd colum 
	os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
	# sort -t$'\t' -k 10,10nr -k 2,2n
	
def sort_single_capture_primers_by_2_columns(In,Out,number_of_primers_to_review, column, column2):
	# this finds the primers that best capture the whole cluster, return the 100 best.
	import os
	os.system("sort -k %i,%inr -k %i,%inr -k 2,2n %s > tempsort.txt" %(column, column, column2, column2,In)) # SORT COMMAND, SORT FIRST ON THE 10th column numerically, next sort on the 2nd colum 
	os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
	# sort -t$'\t' -k 10,10nr -k 2,2n

def sort_by_priority_columns(In,Out,priority_column, second_priority_column, number_of_primers_to_review, direction):
	# this finds the primers that best capture the whole cluster, return the 100 best.
	import os
	priority_column = int(priority_column) 
	second_priority_column = int(second_priority_column)
	if direction is 'F':
		os.system("sort -k %i,%in -k %i,%in %s > tempsort.txt" %(priority_column, priority_column, second_priority_column, second_priority_column, In)) # SORT COMMAND, SORT FIRST ON THE priority column  numerically, next sort on the 2nd priority column 
		os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
	elif direction is 'R':
		os.system("sort -k %i,%inr -k %i,%in %s > tempsort.txt" %(priority_column, priority_column, second_priority_column, second_priority_column, In)) # SORT COMMAND, SORT FIRST ON THE priority column  numerically, next sort on the 2nd priority column 
		os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
	# sort -t$'\t' -k 10,10nr -k 2,2n
	
	
def sort_by_3_priority_columns(In,Out,priority_column, second_priority_column, third_priority_column, number_of_primers_to_review, direction):
	# this finds the primers that best capture the whole cluster, return the 100 best.
	import os
	priority_column = int(priority_column) 
	second_priority_column = int(second_priority_column)
	third_priority_column = int(third_priority_column)
	if direction is 'F':
		os.system("sort -k %i,%in -k %i,%inr -k %i,%in %s > tempsort.txt" %(priority_column, priority_column, second_priority_column, second_priority_column,third_priority_column,third_priority_column, In)) # SORT COMMAND, SORT FIRST ON THE priority column  numerically, next sort on the 2nd priority column 
		os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
	elif direction is 'R':
		os.system("sort -k %i,%inr -k %i,%inr -k %i,%in %s > tempsort.txt" %(priority_column, priority_column, second_priority_column, second_priority_column,third_priority_column,third_priority_column, In)) # SORT COMMAND, SORT FIRST ON THE priority column  numerically, next sort on the 2nd priority column 
		os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
	# sort -t$'\t' -k 10,10nr -k 2,2n


def sort_exclusion_table(In, Out):
	import os
	os.system("sort -k ")
	
	
	
	
	


def probe_plotter(line, fn_blastout, count):
	# ARGS
		# A LINE FROM MY BLAST OUT (NOTE: COLUMN INDEX 12 and -1 MUST CONTAIN THE RELEVANT FIELDS)
		# fn - filename : FOR THE BLASTOUT PUT 
	# DEPDENDS ON grab_sub_blastout and R_probe_plotter.R 
	# RETURNS
		# EXECUTES R, WHICH WRITES A PDF SHOWNING THE PERFORMANCE OF THE ASSAY 
	# RUN ME: probe_plotter(line, '../3_Run/PF13486_Full_Length350to700AA.ncbi.aa.gbk.faa.blastresult')	
		# line = ">BAE84628.1_0	0.054721	GGCTATTATGCAGCGCCGTG	CAGAACTCGCGTACCCCGAA	63.949	64.004	1062	1258	197	1.0	1.0	1.0	BAE84628.1|CAD28792.1	0.005	0.052	0.054	ACH87596.1|CAR57929.1|CAR57926.1|ACH87599.1|ACH87597.1|BAF57046.1|CAR57927.1|AAO60101.1|CAR57937.1|CAR57932.1|CAJ75430.1|CAJ75435.1|CAR57933.1|ACH87594.1|CAR57934.1|ACH87598.1|CAR57936.1|CAR57931.1|CAD28790.2|CAR57935.1|CAR57930.1"
	import os 
	line = line.strip()
	primary_seq = line.split()[0]
	assay_id = primary_seq.replace(">","")
	assay_id = str(count) + "_" + assay_id # THIS ALLOWS COUNTING
	primary_seq = primary_seq.split("_")[0].replace(">","")
	inclusion_hit_list = line.split()[12].split("|")
	exclusion_hit_list = line.split()[16].split("|")
	from pcr_tools import grab_sub_blastout
	grab_sub_blastout(primary_seq, fn_blastout, 'BLASTout.temp') # Grab the sub-section of the blast that you need, next you will add a column with hits and non hits
	fh = open('BLASTout.temp', 'r')
	oh = open('BLASToutScored.temp', 'w')
	for line in fh:
		line = line.strip()
		subject = line.split()[1]
		if subject in inclusion_hit_list:
			oh.write(line + "\tI\n")
		elif subject in exclusion_hit_list:
			oh.write(line + "\tE\n")
		else:
			oh.write(line + "\tN\n")
	fh.close()
	oh.close()
	# CALL R 
	path = os.getcwd() + "/"
	os.system("./R_probe_plotter.R %s %s %s"%(assay_id, path, 'BLASToutScored.temp'))


def advanced_probe_plotter(full_line, fn_blastout, count):
	# ARGS
		# A LINE FROM MY BLAST OUT (NOTE: COLUMN INDEX 12 and -1 MUST CONTAIN THE RELEVANT FIELDS)
		# fn - filename : FOR THE BLASTOUT PUT 
	# DEPDENDS ON grab_sub_blastout and R_complex_probe_plotter.R 
	# RETURNS
		# EXECUTES R, WHICH WRITES A PDF SHOWNING THE PERFORMANCE OF THE ASSAY 
	# RUN ME: 
	
	import os
	full_line = full_line.strip()
	inclusion_hit_list = full_line.split()[10].split("|")
	thermo_penalty = full_line.split()[1] # This is the primer3  penalty 
	inclusivity = full_line.split()[9] # This is the inclusivity of the desired sequences
	primary_seq = full_line.split()[0] # >BAE84628.1_0	
	original_seq = primary_seq.split("_")[0].replace(">","")
	assay_id = primary_seq.replace(">","") #BAE84628.1_0	
	assay_id = str(count) + "_" + assay_id # THIS ALLOWS COUNTING
	primary_seq = primary_seq.split("_")[0].replace(">","") #BAE84628.1
	match_info = full_line.split()[9::2]
	match_hits = full_line.split()[10::2]
	exclusion_hits= match_hits[5:]
	from pcr_tools import list_partition
	non_redundant_hits = list_partition(exclusion_hits, "|", "No_Hits")
	my_dict = dict(enumerate(non_redundant_hits))
	new_dict = {} # Contains the index position of every sequence (in this case 0 corresponds with zero mismatchs, 1 with (1,2 mismatchs), 2 (3,4 mismatches))
	for k in my_dict.keys():
		L = my_dict[k]
		for l in L:
			new_dict[l] = k 
	# NOW WE GO THROUGH THE ACTUAL BLAST AND MARK EACH SEQUENCE WITH LOOKUPS.
	from pcr_tools import grab_sub_blastout
	grab_sub_blastout(primary_seq, fn_blastout, 'BLASTout.temp') # Grab the sub-section of the blast that you need, next you will add a column with hits and non hits
	fh = open('BLASTout.temp', 'r')
	oh = open('BLASToutScored.temp', 'w')
	for line in fh:
	 	line = line.strip()
		subject = line.split()[1]
		if subject in inclusion_hit_list:
			oh.write(line + "\tI" + "\tNA\n")
			continue
		else: 
			oh.write(line + "\tE")
			if subject in new_dict.keys():
				x = new_dict[subject]
				oh.write("\t%s\n" %(x))
			else:
				oh.write("\tNA\n")
	fh.close()
	oh.close()
	path = os.getcwd() + "/"
	
	os.system("../scripts_used/R_complex_probe_plotter.r %s %s %s %s"%(assay_id, path, 'BLASToutScored.temp', original_seq))


def grab_sub_blastout(seq_id, fn_blastout, Out):
	import os 
	os.system("grep -E '^%s' %s > %s" %(seq_id, fn_blastout, Out))


def list_partition(L, deliminator, empty):	
# THIS TAKES A LIST OF STRING WITH A DELIMINATOR BREAKS THEM OPEN AND PARTITIONS INTO A LIST OF LIST, WHERE EACH SUBSEQUENCE OBJECT HAS NO ENTRIES FROM THE PREVIOUS SET. THIS HAS PRACTICAL IMPLICATION WHEN TAKING OUR LIST OF HITS GIVEN INCREASING NUMBER OF ACCEPTABLE MISMATCHES. 
# TEST IT WITH: 
#mega = ["A|B|C","A|B|C|D","A|B|C|D|E","A|B|C|D|E|F"]
#print list_partition(mega, "|", "No_Hits")
# SHOULD GET : [['A', 'C', 'B'], ['D'], ['E'], ['F']]
	def no_hits_check(L, string):
		if string in L :
			return [ ]
		else: 
			return L
	accumulated_set = set()
	output = []
	while len(L) > 0:
		my_list = L.pop(0).split(deliminator)
		my_list = no_hits_check(my_list, empty)
		my_set = set(my_list)
		my_set = my_set - accumulated_set
		accumulated_set = accumulated_set | my_set
		output.append(list(my_set))
	return output







