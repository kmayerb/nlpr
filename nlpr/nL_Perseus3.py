# MODULE FOR USE WITH nL_PRIME.py 
# KEY DATES
# October 10, 2012 
# Tweaked Jan 24, 2012 for use with nL_PRIME.py
def Perseus_part_1(fn_fasta, ref_seq, inclusion_list, number_of_primers_to_generate, number_of_primers_to_review, inclusion_set_threshold):
	print inclusion_list
	from pcr_tools import generate_seq_dictionary
	master_dictionary = generate_seq_dictionary(fn_fasta) # Write a seq dictionary from a fasta, for later use 
	chosen_sequence = str(ref_seq)
	from pcr_tools import write_subfasta
	write_subfasta(inclusion_list,master_dictionary, "inclusion_set.fna") # we pass the List and Dictionary
	#write_subfasta(exclusion_list,master_dictionary, "exclusion_set.fna") # we pass the Second List and Dictionary 
	ref_seq_list = [str(ref_seq)] # make a list of one element, for using write subfasta function
	write_subfasta(ref_seq_list, master_dictionary, "choosen.fna")
	
	####################################
	# NEED TO REMOVE ANY BASES OTHER THAN ATCG (KMB, BUG FIX JAN 30, 2012)
	fh_fix = open("choosen.fna", 'r')
	header = fh_fix.readline()
	line = fh_fix.readline()
	import string
	new_line = line.translate(string.maketrans("RYSWK", "NNNNN")) # GET RID OF MOST COMMON NON_STANDARD AMIGUOUS NUCLEOTIDES
	fh_fix.close()
	fh_fix = fh_fix = open("choosen.fna", 'w')
	fh_fix.write(str(header) + str(new_line))
	fh_fix.close()
	##################################	
	from pcr_tools import fna_to_pr3_v2_3_4_input
	fna_to_pr3_v2_3_4_input("choosen.fna", "choosen.pr3in", number_of_primers_to_generate) # (1) Generate the primer3 input
	from pcr_tools import pr3
	pr3("choosen.pr3in", "choosen.pr3out") # (2) Run Primer 3, produce the primer 3 output file
	from pcr_tools import parse_pr3_v2_3_4_out
	parse_pr3_v2_3_4_out("choosen.pr3out", chosen_sequence + ".pr3set") #(3) generate the output file of primers
	# This block goes through the candidate assays and tests what sequences have complimentary targets
	fh = open(chosen_sequence + ".pr3set", 'r')
	all_scored_primers = open(chosen_sequence + ".scoredprimers" + ".InSet." + str(inclusion_set_threshold), 'w')
	from pcr_tools import hit_percentage_calc_plus_discovery
	for line in fh:
		line = line.strip()
		primer_id = line.split()[0]
		p1 = line.split()[2]
		p2 = line.split()[3]
		from pcr_tools import write_fuznuc_pattern # writes a temporary pattern file.
		write_fuznuc_pattern(p1,p2, 'temp.pat')
		from pcr_tools import execute_fuznuc # passes the pattern file and sequence files to fuznuc. #NOTE: ENSURE THE FOLLOWING IS INSTALLED in ~/EMBOSS-6.5.7/emboss/fuzznuc !!!!!!!!
		execute_fuznuc("temp.pat","inclusion_set.fna", "temp.fzout")
		from pcr_tools import parse_fzout # parse the fuznuc output and populate the dictionary
		D = parse_fzout("temp.fzout")
		#print D 
		X = [1,3,5,6,9] # these integers are the number of mismatches allowed in each round
		output_list = []
		hit_list = []
		all_scored_primers.write(line + "\t")
		for x in X:
			R_Dict = hit_percentage_calc_plus_discovery(inclusion_list,D,x)
			hit_list = R_Dict['hits']
			if float(R_Dict['hit_percentage']) > 0:
				all_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + '|'.join(map(str, hit_list)) + "\t")
			else: 
				all_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + 'No_Hits' + "\t")
		all_scored_primers.write("\n")
	fh.close()
	all_scored_primers.close()
	# This blocks finds good primer pairs that best capture the most sequences in the desired cluster, ranks them by thermodynamic penalty. Note that you might want to instead say anything that captures 90 percent of the clusters will do fine and look for ones with better thermodynamics. 
	from pcr_tools import sort_single_capture_primers_by_2_columns
	sort_single_capture_primers_by_2_columns(chosen_sequence + ".scoredprimers" + ".InSet." + str(inclusion_set_threshold), chosen_sequence + ".best100primers" + ".InSet." + str(inclusion_set_threshold), number_of_primers_to_review, 12, 10)


def Perseus_part_2(fn_fasta, ref_seq, exclusion_list, number_of_primers_to_generate, number_of_primers_to_review, inclusion_set_threshold):
	from pcr_tools import generate_seq_dictionary
	master_dictionary = generate_seq_dictionary(fn_fasta)
	chosen_sequence = str(ref_seq)
	from pcr_tools import write_subfasta
	write_subfasta(exclusion_list, master_dictionary, "exclusion_set.fna")

	fh = open(chosen_sequence + ".best100primers" + ".InSet." + str(inclusion_set_threshold),'r')
	second_scored_primers = open(chosen_sequence + ".scored_exclusion." + "InSet." + str(inclusion_set_threshold) , 'w')
	trigger = 0
	from pcr_tools import hit_percentage_calc_plus_discovery
	while trigger < number_of_primers_to_review - 1:
		trigger = trigger + 1
		line = fh.readline().strip() # the idea here is to read one line at a time as we might not need to go forever. 
		# Get the primer from the line 
		try:
			p1 = line.split()[2]
			p2 = line.split()[3]
			##############KMB FIX Jan 31, 2013. ####################
			# NEED A BREAK COMMAND IF PRIMERS ARE NOT SUFFICIENTLY INCLUSIVE  
			inclusivity_metric_no_MM = line.split()[9]
			inclusivity_metric_two_MM = line.split()[11]
			if float(inclusivity_metric_two_MM) < 0.01 or float(inclusivity_metric_no_MM) < 0.01: # THIS MEANS ALLOWING FOR TWO TOTAL MISMATCHES ALL OF THE INCLUSION SET MUST BE HIT
				continue
			# Write fuznuc pattern
			from pcr_tools import write_fuznuc_pattern
			write_fuznuc_pattern(p1,p2, 'temp.pat', mismatches = 4)
			# Search 
			from pcr_tools import execute_fuznuc
			execute_fuznuc("temp.pat","exclusion_set.fna", "temp.fzout")
			from pcr_tools import parse_fzout
			D = parse_fzout("temp.fzout")
			# import pprint
			# 		pp = pprint.PrettyPrinter(indent=4)
			# 		pp.pprint(D)
			X = [1,3,5,7,9] 
			output_list = []
			hit_list = []
			second_scored_primers.write(line + "\t")
			for x in X:
				R_Dict = hit_percentage_calc_plus_discovery(exclusion_list,D,x) # Needs the list of sequences to exclude
				# output_list.append(R_Dict['hit_percentage'])
				hit_list = R_Dict['hits']
				if float(R_Dict['hit_percentage']) != 0:
					second_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + '|'.join(map(str, hit_list)) + "\t")
				else: 
					second_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + 'No_Hits' + "\t")
			#	second_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + '|'.join(map(str, hit_list))+ "\t") 
				#print hit_list
			second_scored_primers.write("\n")
		except IndexError:
			continue
		#hit_list = list(set(hit_list)) # get uniques only, convert back to a hit_list
		#print hit_list
			#Output the key lines and the hit percentages
		#second_scored_primers.write(line + "\t" + "\t".join(map(str,output_list)) + "\t" + '|'.join(map(str, hit_list)) + "\n")	
	second_scored_primers.close()
	fh.close()


# FOR SITE SPECIFIC CASE, WHERE WE DESIGN FROM A HMM Extracted Fragment, we need to write choosen.fna a little differently 
def Perseus_part_1_hmm_version(fn_fasta, ref_seq, inclusion_list, number_of_primers_to_generate, number_of_primers_to_review, inclusion_set_threshold):
	print inclusion_list
	from pcr_tools import generate_seq_dictionary
	master_dictionary = generate_seq_dictionary(fn_fasta) # Write a seq dictionary from a fasta, for later use 
	chosen_sequence = str(ref_seq)
	from pcr_tools import write_subfasta
	write_subfasta(inclusion_list,master_dictionary, "inclusion_set.fna") # we pass the List and Dictionary
	#write_subfasta(exclusion_list,master_dictionary, "exclusion_set.fna") # we pass the Second List and Dictionary 
	
	###############
	##### MODIFICATION 	
	# REMOVED:
	#ref_seq_list = [str(ref_seq)] # make a list of one element, for using write subfasta function
	#write_subfasta(ref_seq_list, master_dictionary, "choosen.fna")
	# Added: Quickly make choosen.fna from fasta file generated in nL_site_specific_PRIME.py
	import os
	os.system('more temp_fasta_fragment.fna > choosen.fna')
	######################################################
	

	####################################
	# NEED TO REMOVE ANY BASES OTHER THAN ATCG (KMB, BUG FIX JAN 30, 2012)
	fh_fix = open("choosen.fna", 'r')
	header = fh_fix.readline()
	line = fh_fix.readline()
	import string
	new_line = line.translate(string.maketrans("RYSWK", "NNNNN")) # GET RID OF MOST COMMON NON_STANDARD AMIGUOUS NUCLEOTIDES
	fh_fix.close()
	fh_fix = fh_fix = open("choosen.fna", 'w')
	fh_fix.write(str(header) + str(new_line))
	fh_fix.close()
	##################################	
	from pcr_tools import fna_to_pr3_v2_3_4_input
	fna_to_pr3_v2_3_4_input("choosen.fna", "choosen.pr3in", number_of_primers_to_generate) # (1) Generate the primer3 input
	from pcr_tools import pr3
	pr3("choosen.pr3in", "choosen.pr3out") # (2) Run Primer 3, produce the primer 3 output file
	from pcr_tools import parse_pr3_v2_3_4_out
	parse_pr3_v2_3_4_out("choosen.pr3out", chosen_sequence + ".pr3set") #(3) generate the output file of primers
	# This block goes through the candidate assays and tests what sequences have complimentary targets
	fh = open(chosen_sequence + ".pr3set", 'r')
	all_scored_primers = open(chosen_sequence + ".scoredprimers" + ".InSet." + str(inclusion_set_threshold), 'w')
	from pcr_tools import hit_percentage_calc_plus_discovery
	for line in fh:
		line = line.strip()
		primer_id = line.split()[0]
		p1 = line.split()[2]
		p2 = line.split()[3]
		from pcr_tools import write_fuznuc_pattern # writes a temporary pattern file.
		write_fuznuc_pattern(p1,p2, 'temp.pat')
		from pcr_tools import execute_fuznuc # passes the pattern file and sequence files to fuznuc. #NOTE: ENSURE THE FOLLOWING IS INSTALLED in ~/EMBOSS-6.5.7/emboss/fuzznuc !!!!!!!!
		execute_fuznuc("temp.pat","inclusion_set.fna", "temp.fzout")
		from pcr_tools import parse_fzout # parse the fuznuc output and populate the dictionary
		D = parse_fzout("temp.fzout")
		#print D 
		X = [1,3,5,6,9] # these integers are the number of mismatches allowed in each round
		output_list = []
		hit_list = []
		all_scored_primers.write(line + "\t")
		for x in X:
			R_Dict = hit_percentage_calc_plus_discovery(inclusion_list,D,x)
			hit_list = R_Dict['hits']
			if float(R_Dict['hit_percentage']) > 0:
				all_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + '|'.join(map(str, hit_list)) + "\t")
			else: 
				all_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + 'No_Hits' + "\t")
		all_scored_primers.write("\n")
	fh.close()
	all_scored_primers.close()
	# This blocks finds good primer pairs that x capture the most sequences in the desired cluster, ranks them by thermodynamic penalty. Note that you might want to instead say anything that captures 90 percent of the clusters will do fine and look for ones with better thermodynamics. 

### MODIFICATION, I HANDLE SORTING DIRECTLY IN PROGRAM FOR HMM CASE [KMB Jan 11, 2013]
	#from pcr_tools import sort_single_capture_primers
	#sort_single_capture_primers(chosen_sequence + ".scoredprimers" + ".InSet." + str(inclusion_set_threshold), chosen_sequence + ".best100primers" + ".InSet." + str(inclusion_set_threshold), number_of_primers_to_review)


def Perseus_part_2_hmm_version(fn_fasta, ref_seq, exclusion_list, number_of_primers_to_generate, number_of_primers_to_review, inclusion_set_threshold):
	from pcr_tools import generate_seq_dictionary
	master_dictionary = generate_seq_dictionary(fn_fasta)
	
	chosen_sequence = str(ref_seq)
	
	from pcr_tools import write_subfasta
	write_subfasta(exclusion_list, master_dictionary, "exclusion_set.fna")

	fh = open(chosen_sequence + ".best.scoredprimers.InSet." + str(inclusion_set_threshold),'r')
	second_scored_primers = open(chosen_sequence + ".best.scoredprimers.InSet.scored_exclusion." + str(inclusion_set_threshold), 'w')
	
	trigger = 0
	from pcr_tools import hit_percentage_calc_plus_discovery
	while trigger < number_of_primers_to_review - 1:
		trigger = trigger + 1
		line = fh.readline().strip() # the idea here is to read one line at a time as we might not need to go forever. 
		# Get the primer from the line 
		p1 = line.split()[2]
		p2 = line.split()[3]


		##############KMB FIX Jan 31, 2013. ####################
		# NEED A BREAK COMMAND IF PRIMERS ARE NOT SUFFICIENTLY INCLUSIVE  
		inclusivity_metric_no_MM = line.split()[9]
		inclusivity_metric_two_MM = line.split()[11]

		if float(inclusivity_metric_two_MM) < 0.49 or float(inclusivity_metric_no_MM) < 0.49: # THIS MEANS ALLOWING FOR TWO TOTAL MISMATCHES ALL OF THE INCLUSION SET MUST BE HIT
		#### THE INCLUSIVITY METRIC WERE MODIFERED KMB FEb 11, 2013.
		
			continue
		# Write fuznuc pattern
		from pcr_tools import write_fuznuc_pattern
		write_fuznuc_pattern(p1,p2, 'temp.pat', mismatches = 4)
		# Search 
		from pcr_tools import execute_fuznuc
		execute_fuznuc("temp.pat","exclusion_set.fna", "temp.fzout")
		from pcr_tools import parse_fzout
		D = parse_fzout("temp.fzout")
		#import pprint
		#pp = pprint.PrettyPrinter(indent=4)
		#pp.pprint(D)
		X = [1,3,5,7,9] 
		output_list = []
		hit_list = []
		second_scored_primers.write(line + "\t")
		
		for x in X:
			R_Dict = hit_percentage_calc_plus_discovery(exclusion_list,D,x) # Needs the list of sequences to exclude
			# output_list.append(R_Dict['hit_percentage'])
			hit_list = R_Dict['hits']
			if float(R_Dict['hit_percentage']) != 0:
				second_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + '|'.join(map(str, hit_list)) + "\t")
			else: 
				second_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + 'No_Hits' + "\t")
		#	second_scored_primers.write(str(R_Dict['hit_percentage']) + "\t" + '|'.join(map(str, hit_list))+ "\t") 
			#print hit_list
		second_scored_primers.write("\n")
		#hit_list = list(set(hit_list)) # get uniques only, convert back to a hit_list
		#print hit_list
			#Output the key lines and the hit percentages
		#second_scored_primers.write(line + "\t" + "\t".join(map(str,output_list)) + "\t" + '|'.join(map(str, hit_list)) + "\n")	
	second_scored_primers.close()
	fh.close()
