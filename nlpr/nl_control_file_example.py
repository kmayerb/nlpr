#
#	             LL          PP P P PP  RR R R RR  II  MM      MM  EE E E E 
#	             LL          PP     PP  RR     RR  II  MM M  M MM  EE
#	  NN     NN  LL          PP P P PP  RR R R RR  II  MM  M   MM  EE E E
#	  NN N   NN  LL          PP         RR RR      II  MM      MM  EE
#	  NN   N NN  LL          PP         RR   RR    II  MM      MM  EE
#	  NN     NN  LL L L L    PP         RR     RR  II  MM      MM  EE E E E  
#
#
#							Version 1.0 Input File
#							March 2013 (Reused June 2014)
#							Koshlan Mayer-Blackwell
#							Spormann Lab | Stanford University
#
#
# This is a master control file for designing candidate assays for use on a 100nL qPCR platform (With the default parameters, assays should perform @ 60C Ta with the Roche MM on the Wafergen platform) Above the break are editible paths, filenames, and parameters that allow the user to define the program's operations and inputs.
#
#
#
# Paths
#
# Program paths Must Be Set Correctly For Your Computer

# DATA BLOCK 1: PATHS 
	# Blast path 
#blast_path = '/Users/JGrembi/ncbi-blast-2.2.18+/bin'
	# CD-Hit_path 
CDHit_path = '/usr/bin/'#'/usr/bin/cd-hit/'
	# Silix (Requires Boost Installed:  The C++ Boost:program_options package)	)
silix_path = '/software/silix-1.2.6/src/' #'/software/silix-1.2.6/src/silix'
	# Primer3 path
primer3_path = '/usr/bin/' #/usr/bin/primer3_core
	# fuznuc path 
fuznuc_path = '/usr/bin/' #/usr/bin/fuzznuc
	#
utils_path = './utils/' 
temp_path =  './temp/'

# DATA BLOCK 2: GLOBAL VARIABLES AND FILENAMES 
silix_database_name =  'example_SILIX' 
	# primer design run name
run_name = 'example' 
	# filename_fasta_amino_acid
fn_faa ='/nlpr/nlpr/inputs/PF00699.gp.target.faa.minimal.concise' 
	# filename_fasta_nucleic_acid
fn_fna = '/nlpr/nlpr/inputs/PF00699.gp.minimal.nt.gbk.fna.concise' 
	# filename_blast_result
fn_blastresult = '/nlpr/nlpr/inputs/PF00699.gp.target.faa.minimal.concise.all-v-all_blastp_output'  

# custom_ranked_accession_list (ARGUEMNT ONLY USED IF: use_custom_ranked_list ==True)
# The user may provide a custom_ranked list of accession of key proteins to be used as 
# reference for primer design. For instance a known pathogenic isolate can be used as the 
# reference sequence so that primer is guaranteed to match that sequence. If 
# use_custom_ranked_list = False any member of a sequence cluster may be used 
# as the primary reference. 
custom_ranked_accession_list = 	['AAC72748.1']

# Flow Control Arguments 
	# generate_new_silix_database
generate_new_silix_database = True
	# generate_new_subnetwork_run 
generate_new_subnetwork_run = True #True
	# use_custom_ranked_list
use_custom_ranked_list = True
    # use_custom_network_input
custom_network_input = False

	# Generate Primers: DO YOU WANT TO GENERATE PRIMERS
generate_primers_for_A_class = True
generate_primers_for_B_class = True
generate_primers_for_C_class = False
generate_primers_for_D_class = False
generate_primers_for_E_class = False
generate_primers_for_F_class = False


# Default Parameters 
# silix
silix_overlap = 0.8 # Minimum fraction of a sequence making up allignment used to draw a network edge.
# primer3 
number_of_primers_to_generate = 5000	
# fuznuc nl
number_of_primers_to_review = 5000

# Subnetwork parametners
starting_level_identity = 0.90 # 
inclusion_size = 3
inclusion_identity_final_point = 0.87
starting_level_identity2 = 0.90
inclusion_identity_final_point2 = 0.52 # MUST NOT GO LOWER 51
exclusion_size = 10
network_size = 5


######## NETWORK DEFINITION ########
### Assuming you run this from scratch. What is the computer doing? First you must generate a silix database which is primarily done by running the program Silix repeatedly on a blastp output file with various threshold values
import os 
import sys
from nL_network_tools import print_nL_prime_logo

if generate_new_silix_database == True:
	from nL_network_tools import create_silix_database
	create_silix_database(silix_database_name , silix_path, silix_overlap, fn_faa, fn_blastresult)
	sys.stderr.write("\nYou wrote the -- %s.silix.archive --\n\tusing <nL_network_tools.create_silix_database>\n" %(silix_database_name))
print_nL_prime_logo()
if generate_new_subnetwork_run ==True:
	from nL_network_tools import establish_rankings
	if use_custom_ranked_list == False:
		ranked_accession_list = establish_rankings(silix_database_name, 0.97)
		sys.stderr.write("\nYou ranked accessions in -- %s.silix.archive --\n\tand passed the ranked list to <nL_network_tools.determine_sub_networks>\n" %(silix_database_name))
	else:
		ranked_accession_list = custom_ranked_accession_list
	
	from nL_network_tools import determine_sub_networks
	determine_sub_networks(silix_database_name ,run_name, ranked_accession_list,inclusion_size,starting_level_identity,starting_level_identity2, inclusion_identity_final_point,inclusion_identity_final_point2, network_size, exclusion_size)
	sys.stderr.write("\nYou wrote -- %s.sub_networks --\n\tusing <nL_network_tools.determine_sub_networks>\n\n" %(silix_database_name))

	from nL_network_tools import summarize_sub_network_performance
	summarize_sub_network_performance(silix_database_name, run_name)
	import time
	time.sleep(2)


##### PRIMER DESIGN AND REPORTING	
if generate_primers_for_A_class == True:
	from nL_network_tools import sub_file_by_grade_x
	grade = 'A'
	sub_file_by_grade_x(run_name,grade) # run_name.sub_networks => run_name.sub_networks.temp
	fn_subnetwork_file = run_name + '.sub_networks.temp'
	
	if not os.path.isdir(run_name):
		os.system('mkdir %s' %(run_name))
	final_output_path = './%s/'%(run_name)
	from nL_network_tools import nL_PRIME_A
	
	nL_PRIME_A(fn_fna, fn_blastresult, fn_subnetwork_file, number_of_primers_to_generate, number_of_primers_to_review, utils_path, temp_path, final_output_path, grade)

if generate_primers_for_B_class == True:
	from nL_network_tools import sub_file_by_grade_x
	grade = 'B'
	sub_file_by_grade_x(run_name,grade) # run_name.sub_networks => run_name.sub_networks.temp
	fn_subnetwork_file = run_name + '.sub_networks.temp'
	if not os.path.isdir(run_name):
		os.system('mkdir %s' %(run_name))
	final_output_path = './%s/'%(run_name)
	from nL_network_tools import nL_PRIME_B
	grade = 'B'
	nL_PRIME_B(fn_fna, fn_blastresult, fn_subnetwork_file, number_of_primers_to_generate, number_of_primers_to_review, utils_path, temp_path, final_output_path, grade)
	
if generate_primers_for_C_class == True:
	pass
if generate_primers_for_D_class == True:
	pass
if generate_primers_for_E_class == True:
	pass
if generate_primers_for_F_class == True:
	pass 
#### TOMORROW YOU CAN ADD nL_prime_site specific functionality to this.




