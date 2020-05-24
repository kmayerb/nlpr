# Koshlan Mayer-Blackwwell 
# March 8, 2013
# The following tools are used to determine clusters/sub_networks of related proteins.

def run_silix(silix_path, identity, overlap, fn_fasta, fn_blastresult):
	'''
	Method name: run_silix
	march 8 2013 *kmb
	
	Description: Executes the Clustering Program Silix from the University of Lyon, France. http://lbbe.univ-lyon1.fr/Documentation,3012.html?lang=en 
	
	Inputs:
		--silix_path--			path to program
	 	--identity -- 			float "Min  % identity to accept blast hits for building families (in [0,1], 0.35 by default )"
		--overlap--				float Min % overlap to accept blast hits for  building  families (in [0,1], 0.80 by default)
		--fn_fasta-- 			fasta filename
		--fn_blastresult--		reciprocal tabular blast from file
		
		
	Writes:
		--silix_out.temp --	file with accesion and cluster number
	Stores in object attributes:
		--  --
	'''
	import os
	os.system('%ssilix -i %f -r %f %s %s > silix_out.temp' %(silix_path, identity, overlap, fn_fasta, fn_blastresult))


def parse_silix(fn_in = 'silix_out.temp'):
	'''
	Method name: pares_silix
	march 8 2013 *kmb
	
	Description: Silix produces a simple output of one sequences accession and cluster_number per line. 
	This method parses that informatoin into two related dictionaries
			-- dict_cluster_to_accession_list --	keys are cluster to a list of accessoins D[cluster_1] = [seq1,seq2,seq3]
			-- dict_accession_to_cluster -- 		keys are accesions  D[seq1] = cluster_1
	These dictionaries are useful because you can use them to link one sequence to all of the sequences in the cluster
	
	Inputs:
		-- fn_in -- 	 the filename of the silex output files to be parsed 
	Returns: 
		-- dict_accession_to_cluster --			keys are accesions  e.g. D[seq1] = cluster_1
		-- dict_cluster_to_accession_list --	keys are cluster to a list of accessoins e.g. D[cluster_1] = [seq1,seq2,seq3]
	'''
	dict_accession_to_cluster = dict()
	dict_cluster_to_accession_list = dict()
	fh = open(fn_in,'r')
	for line in fh:
		cluster, accession = line.strip().split("\t")
		dict_accession_to_cluster[accession] = cluster
		if cluster not in dict_cluster_to_accession_list.keys():
			dict_cluster_to_accession_list[cluster] = list()
			dict_cluster_to_accession_list[cluster].append(accession)
		else: 
			dict_cluster_to_accession_list[cluster].append(accession)
	return dict_accession_to_cluster, dict_cluster_to_accession_list


def create_silix_database(database_name, silix_path, overlap, fn_fasta, fn_blastresult):
	'''
	Method name: create_silix_database
	march 8 2013 *kmb
	
	Description: When running Silix for many possibe edge threholds (sequence idenity from 0.99 - 0.50) with default overlap (.50) 
	This method generates a Master Dictionary allowing one to find how many nodes are connected to any node of interest at any level of similiarity
	
	Inputs:
		-- database_name -- 	this is need to name and store the pickled silix databse
		-- silix_path-- 		this must be passed to < run_silix >
		-- overlap -- 			this must be passed to < run_silix > 
		-- fn_fasta --  		this must be passed to < run_silix > 
		-- fn_blastresult -- 	this must be passed to < run_silix >
	Writes:
		-- 
		-- database_name.silix.archive --			pickled dictionary database dictionary of list of dictionary MD[identiy]= list([accession_to_cluster][cluster_to_list_of_neighbors])
	'''
	master_dictionary = dict()
	identities=[0.99,0.98,0.97,0.96,0.95,0.94,0.93,0.92,0.91,0.9,0.89,0.88,0.87,0.86,0.85,0.84,0.83,0.82,0.81,0.8,0.79,0.78,0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.7,0.69,0.68,0.67,0.66,0.65,0.64,0.63,0.62,0.61,0.6,0.59,0.58,0.57,0.56,0.55,0.54,0.53,0.52,0.51]
	for i in identities:
		run_silix(silix_path, i, overlap, fn_fasta, fn_blastresult)
		D1, D2 = parse_silix()
		master_dictionary[i] = list()
		master_dictionary[i].append(D1)
		master_dictionary[i].append(D2)
		
		import pickle
		fh = open(database_name + ".silix.archive",'w')
		pickle.dump(master_dictionary, fh)
		fh.close()
		


def establish_rankings(database_name,id_level):
	''' given a idenity percentage return a ranked list of the sequences form most to least connected'''
	import pickle
	pickleFileHandle = open(database_name+ ".silix.archive")
	master_dictionary = pickle.load(pickleFileHandle)
	pickleFileHandle.close()
	unranked_list = list()
	Da_to_c, Dc_to_a = master_dictionary[id_level]
	
	for acc in Da_to_c.keys():
		t = acc, len(Dc_to_a[Da_to_c[acc]])	
		unranked_list.append(t)		
	ranked_list = sorted(unranked_list, key=lambda tup: tup[1], reverse = True)
	#for x in ranked_list:
	#	print x[0]
	return [x[0] for x in ranked_list] # kicks out a ranked list of accesions
	


def HS_determine_sub_networks(database_name, run_name, list_rankings, id_level):
	'''
	This was modified from the method < determine_sub_networks() > especially for the high sensitivity (HS) version of the program. 
	The major changes are:
	(1) The addition of a id_level argument that tells the algorithm at what edge threshold to start (0.99 instead of 0.97 in the regular version)
	(2) Whereas the prevoius version of the program kept droping the edge threshold until it found a network of size three or reached a lower edge threshold condition,
	the HS version does not drop below < id_level >
	
	Inputs:
		-- database_name -- 	this is what is used for accessing the silix.archive database and for naming output
		-- list_rankings -- 	this is what is used for poping reference sequences
	
	Writes:
		--database_name.sub_networks -- 
	
	'''	
	def accept_drop_or_fail(ref_seq, inclusion_identity, network_size, acceptable_size = 3):
		''' A sub routine that searches for the size of the network linked to a < ref_seq > given an <inclusion_identity> threshold
		it runs within the context of sub_networks where a < master_dictionary > has already been loaded, see above desciptions for structure of the dual dictionaries
		If the network size is less than the <acceptable_size> the function goes recursive, droping the inclusion identity threshold'''
		Da_to_c, Dc_to_a = master_dictionary[inclusion_identity] # GET THE DICTIONARIES ASSOCIATED WITH THAT IDENTITY LEVEL
		networksize = len(Dc_to_a[Da_to_c[ref_seq]]) # ASK FOR THE NETWORK SIZE, OF
		if inclusion_identity < id_level: # If the inclusion identity drops to low
			return 'failed_inclusion', networksize, inclusion_identity
		if networksize >= acceptable_size: # 
			return 'passed_inclusion', networksize, inclusion_identity
		else:
			new_inclusion_identity = inclusion_identity - 0.01 # Drop the identity and try agains
			new_inclusion_identity = round(new_inclusion_identity, 2)
			return accept_drop_or_fail(ref_seq, new_inclusion_identity, networksize, acceptable_size = 3)
	
	def accept_drop_or_fail2(ref_seq, inclusion_identity, network_size, acceptable_size = 3):
		''' A second sub routine similar to above accept_drop_or_fail, but witha different minimum inclusion identity. 
		It that searches for the size of the network linked to a < ref_seq > given an <inclusion_identity> threshold
		it runs within the context of sub_networks where a < master_dictionary > has already been loaded, see above desciptions for structure of the dual dictionaries
		If the network size is less than the <acceptable_size> the function goes recursive, droping the inclusion identity threshold'''
		
		Da_to_c, Dc_to_a = master_dictionary[inclusion_identity] # GET THE DICTIONARIES ASSOCIATED WITH THAT IDENTITY LEVEL
		networksize = len(Dc_to_a[Da_to_c[ref_seq]]) # ASK FOR THE NETWORK SIZE
		if inclusion_identity < 0.52: # If the inclusion identity drops to low
			return 'failed_exclusion', networksize, inclusion_identity
		if networksize > acceptable_size: # 
			return 'passed_exclusion', networksize, inclusion_identity
		else:
	
			new_inclusion_identity = inclusion_identity - 0.01 # Drop the identity and try against
			new_inclusion_identity = round(new_inclusion_identity, 2)
			return accept_drop_or_fail2(ref_seq, new_inclusion_identity, networksize, acceptable_size)
	import pickle
	pickleFileHandle = open(database_name + ".silix.archive")
	master_dictionary = pickle.load(pickleFileHandle)
	pickleFileHandle.close()
	output_handle = open(run_name + ".sub_networks", 'w')
	while list_rankings:
		starting_level_identity = id_level # Choose at this level to avoid collecting simply near duplicate sequences
		starting_network_size = 0
		ref_seq = list_rankings.pop(0) # GET THE HIGHEST RANKED SEQUENCE 
		result, nsize, inclusion_identity = accept_drop_or_fail(ref_seq, starting_level_identity, starting_network_size, 3)	
		# GET THE INCLUDED SEQUENCES AND REMOVE THEM FROM THE FULL LIST
		Da_to_c, Dc_to_a = master_dictionary[inclusion_identity] # GET THE DICTIONARIES ASSOCIATED WITH THAT IDENTITY LEVEL
		inclusion_set = Dc_to_a[Da_to_c[ref_seq]]
	
		inclusion_set_size = len(inclusion_set)
		starting_level_identity = 0.70
		starting_network_size = 0
	
		result2, nsize2, inclusion_identity2 = accept_drop_or_fail2(ref_seq, starting_level_identity, starting_network_size, inclusion_set_size)	
		Da_to_c, Dc_to_a = master_dictionary[inclusion_identity2] # GET THE DICTIONARIES ASSOCIATED WITH THAT IDENTITY LEVEL
		inclusion_set2 = Dc_to_a[Da_to_c[ref_seq]]
	
		exclusion_set = list(set(inclusion_set2) - set(inclusion_set))
		# OUTPUT WRITE
		# EXPECTED FLAT FILE from < nL_sub_network_heuristic.py >
			# AAW40342.1	8	5	3	94	50	ABA64547.1|AAW40342.1|ADC74667.1|ACZ62492.1|CAI83612.1	AAW40361.1|CAI83519.1|ABY28310.1
			# Column 1 		Accession
			# Column 2		# Full Network Set Size
			# Column 3		# Specific Inclusion Set Size
			# Column 4		# Specific Exclusion Set Size 
			# Column 5		% similarity edge threshold used to get > 2 sequences inclusion set
			# Column 6		% similarity edge threshold used to get atleast 1 sequence in the exclusion
			# Column 7 		Accessions of Inclusion Set "|" -DELIMITED
			# Column 8		Accessions of Exclusion Set "|" -DELIMITED
		output_handle.write(
		ref_seq +"\t"+ 							# Column 1 		Accession
		str(nsize2)+"\t"+						# Column 2		Full Network Set Size
		str(nsize)+"\t"+						# Column 3		Specific Inclusion Set Size
		str(nsize2-nsize)+"\t"+					# Column 4		# Specific Exclusion Set Size 
		str(inclusion_identity)+"\t"+			# Column 5		% similarity edge threshold used to get > 2 sequences inclusion set
		str(inclusion_identity2)+"\t"+			# Column 6		% similarity edge threshold used to get atleast 1 sequence in the exclusion
		"|".join(map(str,inclusion_set))+"\t"+	# Column 7 		Accessions of Inclusion Set "|" -DELIMITED
		"|".join(map(str,exclusion_set))+"\n")	# Column 8		Accessions of Exclusion Set "|" -DELIMITED
	
		# REMOVE THE MEMBERS OF THE INCLUSION SET FROM CONSIDERATION AS POTENTIAL REFERENCE SEQUENCES	
		for x in inclusion_set:
			try:
				list_rankings.remove(x) #REMOVE THE INCLUSION SET FROM THE FULL LIST
			except ValueError:
				pass
	output_handle.close()
	


def determine_sub_networks(database_name, run_name, list_rankings, inclusion_size, starting_level_identity, starting_level_identity2 , inclusion_identity_final_point ,inclusion_identity_final_point2, network_size, exclusion_size):
	'''
	Inputs:
		-- database_name -- 	this is what is used for accessing the silix.archive database and for naming output
		-- list_rankings -- 	this is what is used for poping reference sequences
	
	Writes:
		--database_name.sub_networks -- 
	'''	
	def accept_drop_or_fail(ref_seq, inclusion_identity, acceptable_size, inclusion_identity_final_point):
		''' A sub routine that searches for the size of the network linked to a < ref_seq > given an <inclusion_identity> threshold
		it runs within the context of sub_networks where a < master_dictionary > has already been loaded, see above desciptions for structure of the dual dictionaries
		If the network size is less than the <acceptable_size> the function goes recursive, droping the inclusion identity threshold
		    # Ref_seq
		    # inclusion_identity
		    # network_size
		    # acceptable_size
		    # inclusion_identity_final_point
		'''
		# Say inclusion_identity is 0.97, Acces the dictionary entry at that level
		# See how big the network size is for ref_seq at 0.97
		# Then you check two criteria (have we dropped too low) and (have we found an acceptably large network)
		# If we we haven't then we must recursively try again at a lower inclusion identity
		print acceptable_size,
		Da_to_c, Dc_to_a = master_dictionary[inclusion_identity] # GET THE DICTIONARIES ASSOCIATED WITH THAT IDENTITY LEVEL
		networksize = len(Dc_to_a[Da_to_c[ref_seq]]) # ASK FOR THE NETWORK SIZE, OF
		if inclusion_identity < inclusion_identity_final_point: # If the inclusion identity drops to low
			return 'failed_inclusion', networksize, inclusion_identity
		if networksize >= acceptable_size: # 
			return 'passed_inclusion', networksize, inclusion_identity
		else:
			# RECURSIVELY TRY AGAIN
			new_inclusion_identity = inclusion_identity - 0.01 # Drop the identity and try agains
			new_inclusion_identity = round(new_inclusion_identity, 2)
			return accept_drop_or_fail(ref_seq, new_inclusion_identity, acceptable_size, inclusion_identity_final_point)
	
	
	
	import pickle
	pickleFileHandle = open(database_name + ".silix.archive")
	master_dictionary = pickle.load(pickleFileHandle)
	pickleFileHandle.close()
	output_handle = open(run_name + ".sub_networks", 'w')
	
	
	while list_rankings:
		my_ref_seq = list_rankings.pop(0) # GET THE HIGHEST RANKED SEQUENCE 
		my_starting_level_identity = starting_level_identity # Choose at this level to avoid collecting simply near duplicate sequences
		my_acceptable_size = inclusion_size
		my_inclusion_identity_final_point = inclusion_identity_final_point 
		##### ! ######
		# Ref_seq
		# inclusion_identity
		# acceptable_size
		# inclusion_identity_final_point
		result, nsize, inclusion_identity = accept_drop_or_fail(my_ref_seq, my_starting_level_identity, my_acceptable_size, my_inclusion_identity_final_point)	
		# GET THE INCLUDED SEQUENCES AND REMOVE THEM FROM THE FULL LIST
		Da_to_c, Dc_to_a = master_dictionary[inclusion_identity] # GET THE DICTIONARIES ASSOCIATED WITH THAT IDENTITY LEVEL
		inclusion_set = Dc_to_a[Da_to_c[my_ref_seq]]
		my_starting_level_identity = starting_level_identity2
		my_inclusion_set_size = len(inclusion_set)
		my_acceptable_size = my_inclusion_set_size + exclusion_size
		my_inclusion_identity_final_point = inclusion_identity_final_point2
		##### ! ######	
		result2, nsize2, inclusion_identity2 = accept_drop_or_fail(my_ref_seq, my_starting_level_identity, my_acceptable_size, my_inclusion_identity_final_point)
		Da_to_c, Dc_to_a = master_dictionary[inclusion_identity2] # GET THE DICTIONARIES ASSOCIATED WITH THAT IDENTITY LEVEL
		inclusion_set2 = Dc_to_a[Da_to_c[my_ref_seq]]
		
		exclusion_set = list(set(inclusion_set2) - set(inclusion_set))
		# OUTPUT WRITE
		# EXPECTED FLAT FILE from < nL_sub_network_heuristic.py >
			# AAW40342.1	8	5	3	94	50	ABA64547.1|AAW40342.1|ADC74667.1|ACZ62492.1|CAI83612.1	AAW40361.1|CAI83519.1|ABY28310.1
			# Column 1 		Accession
			# Column 2		# Full Network Set Size
			# Column 3		# Specific Inclusion Set Size
			# Column 4		# Specific Exclusion Set Size 
			# Column 5		% similarity edge threshold used to get > 2 sequences inclusion set
			# Column 6		% similarity edge threshold used to get atleast 1 sequence in the exclusion
			# Column 7 		Accessions of Inclusion Set "|" -DELIMITED
			# Column 8		Accessions of Exclusion Set "|" -DELIMITED
		output_handle.write(
		my_ref_seq +"\t"+ 							# Column 1 		Accession
		str(nsize2)+"\t"+						# Column 2		Full Network Set Size
		str(nsize)+"\t"+						# Column 3		Specific Inclusion Set Size
		str(nsize2-nsize)+"\t"+					# Column 4		# Specific Exclusion Set Size 
		str(inclusion_identity)+"\t"+			# Column 5		% similarity edge threshold used to get > 2 sequences inclusion set
		str(inclusion_identity2)+"\t"+			# Column 6		% similarity edge threshold used to get atleast 1 sequence in the exclusion
		"|".join(map(str,inclusion_set))+"\t"+	# Column 7 		Accessions of Inclusion Set "|" -DELIMITED
		"|".join(map(str,exclusion_set))+"\n")	# Column 8		Accessions of Exclusion Set "|" -DELIMITED
				
		# REMOVE THE MEMBERS OF THE INCLUSION SET FROM CONSIDERATION AS POTENTIAL REFERENCE SEQUENCES	
		for x in inclusion_set:
			try:
				list_rankings.remove(x) #REMOVE THE INCLUSION SET FROM THE FULL LIST
			except ValueError:
				pass
	output_handle.close()


def cytoscape_visualize(fn_in):
	oh_all_reference = open(fn_in + ".cytolist_all_ref", 'w')
	oh_all_inclusion = open(fn_in + ".cytolist_all_inclusion", 'w')
	oh_all = open(fn_in + ".cytolist_all", 'w')
	fh = open(fn_in, 'r')
	for line in fh:
		line = line.strip().split("\t")
		inclusion_number = int(line[1])
		broad_number = int(line[2])
		exclusion_number = int(broad_number) - int(inclusion_number)
		if inclusion_number >= 3 and exclusion_number >= 1:
			ref_seq = line[0]
			inclusion_list =line[-2].split("|")
			broad_list =line[-1].split("|")
			exclusion_list = list(set(broad_list) - set(inclusion_list))
			oh_all_reference.write('%s\n'%(ref_seq))
			[oh_all_inclusion.write(x+"\n") for x in inclusion_list]
		
			oh_all.write('%s # ref_seq\n\n'%(ref_seq))
			[oh_all.write(x + " # inclusion_set_" + ref_seq + "\n") for x in inclusion_list]
			oh_all.write("\n")
		  	[oh_all.write(x + " # exclusion_set_" + ref_seq + "\n") for x in exclusion_list]
			oh_all.write("\n")
	fh.close()		
	oh_all_reference.close()
	oh_all_inclusion.close()
	oh_all.close()	

	
def summarize_sub_network_performance(database_name, run_name):
	fh =open(run_name+'.sub_networks', 'r')
	oh =open(run_name +'.sub_networks.summary', 'w')
	holder_dictionary = { 'A':[list(),list()], 'B':[list(),list()], 'C':[list(),list()], 'D':[list(),list()], 'E':[list(),list()], 'F':[list(),list()] }
	
	for line in fh:
		line = line.strip().split("\t")
		ref_seq = str(line[0])
		inclusion_number = int(line[2])
		broad_number = int(line[1])
		exclusion_number = int(broad_number) - int(inclusion_number)
		inclusion_exclusion_thresholds = (float(line[4]), float(line[5]))
		if inclusion_number >= 3 and exclusion_number >= 1:
			holder_dictionary['A'][0].append(ref_seq)
			holder_dictionary['A'][1].append(inclusion_exclusion_thresholds)
		elif inclusion_number >= 3 and exclusion_number == 0 :
			holder_dictionary['B'][0].append(ref_seq)
			holder_dictionary['B'][1].append(inclusion_exclusion_thresholds)
		elif inclusion_number == 2 and exclusion_number >= 1:
			holder_dictionary['C'][0].append(ref_seq)
			holder_dictionary['C'][1].append(inclusion_exclusion_thresholds)
		elif inclusion_number == 2 and exclusion_number == 0:
			holder_dictionary['D'][0].append(ref_seq)
			holder_dictionary['D'][1].append(inclusion_exclusion_thresholds)
		elif inclusion_number == 1 and exclusion_number >= 1:
			holder_dictionary['E'][0].append(ref_seq)
			holder_dictionary['E'][1].append(inclusion_exclusion_thresholds)
		elif inclusion_number == 1 and exclusion_number == 0:
			holder_dictionary['F'][0].append(ref_seq)
			holder_dictionary['F'][1].append(inclusion_exclusion_thresholds)
	fh.close()
	
	print "\nSummary of Results:\n"
	print "The number of A(+/+)		clusters inclusion_number >= 3 and exclusion_number >= 1"
	print "\t %i\n" %(len(holder_dictionary['A'][0]))
	print "The number of B(+/-) 		clusters inclusion_number >= 3 and exclusion_number == 0"
	print "\t %i\n" %(len(holder_dictionary['B'][0]))
	print "The number of C(+/+) 		clusters inclusion_number == 2 and exclusion_number >= 1"
	print "\t %i\n" %(len(holder_dictionary['C'][0]))
	print "The number of D(+/-) 		clusters inclusion_number == 2 and exclusion_number == 0"
	print "\t %i\n" %(len(holder_dictionary['D'][0]))
	print "The number of E(-/+) 		singles inclusion_number == 1 and exclusion_number >= 1"
	print "\t %i\n" %(len(holder_dictionary['E'][0]))
	print "The number of F(-/-) 		singles inclusion_number == 1 and exclusion_number == 0"
	print "\t %i\n" %(len(holder_dictionary['F'][0]))
	print "\nSummary of Results:\n"
	
	print "Ref_Seq for  A(+/+)			clusters inclusion_number >= 3 and exclusion_number >= 1"
	for x in range(len(holder_dictionary['A'][0])):
		print holder_dictionary['A'][0][x] + "\t" + "|".join(map(str, holder_dictionary['A'][1][x]))
	print "Ref_Seq for  B(+/-) 		clusters inclusion_number >= 3 and exclusion_number == 0"
	for x in range(len(holder_dictionary['B'][0])):
		print holder_dictionary['B'][0][x] + "\t" + "|".join(map(str, holder_dictionary['B'][1][x]))
	print "Ref_Seq for  C(+/+) 		clusters inclusion_number == 2 and exclusion_number >= 1"
	for x in range(len(holder_dictionary['C'][0])):
		print holder_dictionary['C'][0][x] + "\t" + "|".join(map(str, holder_dictionary['C'][1][x]))
	print "Ref_Seq for  D(+/-) 		clusters inclusion_number == 2 and exclusion_number == 0"
	for x in range(len(holder_dictionary['D'][0])):
		print holder_dictionary['D'][0][x] + "\t" + "|".join(map(str, holder_dictionary['D'][1][x]))
	print "Ref_Seq for  E(-/+) 		singles inclusion_number == 1 and exclusion_number >= 1"
	for x in range(len(holder_dictionary['E'][0])):
		print holder_dictionary['E'][0][x] + "\t" + "|".join(map(str, holder_dictionary['E'][1][x]))
	print "Ref_Seq for  F(-/-) 			singles inclusion_number == 1 and exclusion_number == 0"
	for x in range(len(holder_dictionary['F'][0])):
		print holder_dictionary['F'][0][x] + "\t" + "|".join(map(str, holder_dictionary['F'][1][x]))
	oh.write( "\nSummary of Results:\n")
	oh.write( "The number of A(+/+)		clusters inclusion_number >= 3 and exclusion_number >= 1")
	oh.write( "\t %i\n" %(len(holder_dictionary['A'][0])) )
	oh.write( "The number of B(+/-) 		clusters inclusion_number >= 3 and exclusion_number == 0")
	oh.write( "\t %i\n" %(len(holder_dictionary['B'][0])) )
	oh.write( "The number of C(+/+) 		clusters inclusion_number == 2 and exclusion_number >= 1")
	oh.write( "\t %i\n" %(len(holder_dictionary['C'][0])) )
	oh.write( "The number of D(+/-) 		clusters inclusion_number == 2 and exclusion_number == 0")
	oh.write( "\t %i\n" %(len(holder_dictionary['D'][0])) )
	oh.write( "The number of E(-/+) 		singles inclusion_number == 1 and exclusion_number >= 1")
	oh.write( "\t %i\n" %(len(holder_dictionary['E'][0])) )
	oh.write( "The number of F(-/-) 		singles inclusion_number == 1 and exclusion_number == 0")
	oh.write( "\t %i\n" %(len(holder_dictionary['F'][0])) )
	oh.write( "\nSummary of Results:\n")
	
	oh.write( "Ref_Seq for  A(+/+)			clusters inclusion_number >= 3 and exclusion_number >= 1")
	for x in range(len(holder_dictionary['A'][0])):
		oh.write( holder_dictionary['A'][0][x] + "\t" + "|".join(map(str, holder_dictionary['A'][1][x])) )
	oh.write( "Ref_Seq for  B(+/-) 		clusters inclusion_number >= 3 and exclusion_number == 0")
	for x in range(len(holder_dictionary['B'][0])):
		oh.write( holder_dictionary['B'][0][x] + "\t" + "|".join(map(str, holder_dictionary['B'][1][x])) )
	oh.write( "Ref_Seq for  C(+/+) 		clusters inclusion_number == 2 and exclusion_number >= 1")
	for x in range(len(holder_dictionary['C'][0])):
		oh.write( holder_dictionary['C'][0][x] + "\t" + "|".join(map(str, holder_dictionary['C'][1][x])) )
	oh.write( "Ref_Seq for  D(+/-) 		clusters inclusion_number == 2 and exclusion_number == 0")
	for x in range(len(holder_dictionary['D'][0])):
		oh.write( holder_dictionary['D'][0][x] + "\t" + "|".join(map(str, holder_dictionary['D'][1][x])) )
	oh.write( "Ref_Seq for  E(-/+) 		singles inclusion_number == 1 and exclusion_number >= 1")
	for x in range(len(holder_dictionary['E'][0])):
		oh.write( holder_dictionary['E'][0][x] + "\t" + "|".join(map(str, holder_dictionary['E'][1][x])) )
	oh.write( "Ref_Seq for  F(-/-) 			singles inclusion_number == 1 and exclusion_number == 0")
	for x in range(len(holder_dictionary['F'][0])):
		oh.write( holder_dictionary['F'][0][x] + "\t" + "|".join(map(str, holder_dictionary['F'][1][x])) )
	oh.close()
	


def sub_file_by_grade_x(run_name, grade):
	'''
	Take a run_file and only accept those that meet a certain <grade> [A-B-C-D-E-F] make a temporary outfile
	'''
	fh = open(run_name + '.sub_networks','r') # FULL FILE
	oh = open(run_name + '.sub_networks.temp','w') # SUBFILE BY GRADE A
	
	def assessment(inclusion_number, exclusion_number):
		''' Grades information for primer design'''
		exclusion_number = int(exclusion_number)
		inclusion_number = int(inclusion_number)
		if inclusion_number >= 3 and exclusion_number >= 1:
			return 'A'
		elif inclusion_number >= 3 and exclusion_number == 0 :
			return 'B'
		elif inclusion_number == 2 and exclusion_number >= 1:
			return 'C'
		elif inclusion_number == 2 and exclusion_number == 0:
			return 'D'
		elif inclusion_number == 1 and exclusion_number >= 1:
			return 'E'
		elif inclusion_number == 1 and exclusion_number == 0:
			return 'F'
	
	for line in fh:
		print line
		try: 
			ref_seq, broad_count, inclusion_count, exclusion_count, include_threshold, broad_threshold, inclusion_list, exclusion_list = line.strip().split("\t")
		except ValueError: # for those cases with no exclusion set
			ref_seq, broad_count, inclusion_count, exclusion_count, include_threshold, broad_threshold, inclusion_list = line.strip().split("\t")
		if assessment(inclusion_count, exclusion_count) == grade:
			oh.write(line)
		else:
			pass
		
	fh.close()
	oh.close()


def move_all_files_starting_with_x_to_destination(destination_path, x, mypath = "./"):
    '''To keep things neat and tidy this function was written to move all files starting with a string from one path to another
    '''
    import os
    import os.path
    onlyfiles = [ f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath,f)) ]
    onlyfiles_matching = [i for i in onlyfiles if i.startswith(x)]
    for i in onlyfiles_matching:
        os.system('mv %s %s%s' %(i,destination_path,i))


def nL_PRIME_A(fn_nucleotide_fasta, fn_blastresult, fn_subnetwork_file, number_of_primers_to_generate, number_of_primers_to_review,utils_path, temp_path, final_output_path,  grade): # [probably need to feed more paths down the layers of the program]
	'''
	Name: nL_PRIME_A
	
	Description: This is the primer design pipeline for those A-type reference sequences that have >= 3 inclusion set sequences and atleast on exclusion set sequence.
	
	Inputs:
		-- fn_nucleotide_fasta -- 
		-- fn_blastresult --
		-- fn_subnetwork_file --
		-- final_output_path --
		-- number_of_primers_to_generate -- 
		-- number_of_primers_to_review -- 
		-- final_output_path -- 
		-- utils_path --
		-- temp_path --  
		-- grade -- 
		
	Results :
		-- -- Major results is the primers to order file generated in the run_name labeled folder, which contains specific and extended primers that performed best and worst against the penalty function. 
	'''
	import sys
	import os 
	import datetime
	now = datetime.datetime.now()
	appended_final_output_name = 'nL_PRIMER_%s_output_for_orders_'%(grade)+ now.strftime("%Y-%m-%d_%H:%M")
	fh = open(fn_subnetwork_file, 'r')
	for line in fh:
		line = line.strip()
		units = line.split('\t')
		reference_sequence = units[0]
		print "#### REFERENCE : " +  str(reference_sequence)
		specific_exclusion_set = units[-1].split("|")
		if len(specific_exclusion_set) == 0:
			specific_exclusion_set = ['ACZ62391.1'] # DUMMY EXCLUSION SET
		specific_inclusion_set = units[-2].split("|")
		inclusion_set_threshold = units[-4]
		## DETERMINATION OF THE CANDIDATE ASSAYS THAT UNIVERSALLY MATCH THE INCLUSION SET
		from nL_Perseus3 import Perseus_part_1
		Perseus_part_1(fn_nucleotide_fasta, reference_sequence, specific_inclusion_set, number_of_primers_to_generate, number_of_primers_to_review, inclusion_set_threshold)
		##### NOTE: THE FLAT FILE PRODUCED BY THIS BLOCK IS def IS (self, : ):
			# <ref_seq>.BAE84628.1.scoredprimers.InSet.<inclusion_set_threshold> 
			# <ref_seq>.BAE84628.1.BAE84628.1.best100primers.InSet.<inclusion_set_threshold> 
		## SCORE THOSE ASSAYS AGAINST THE EXCLUSION SET
		from nL_Perseus3 import Perseus_part_2
		Perseus_part_2(fn_nucleotide_fasta, reference_sequence, specific_exclusion_set, number_of_primers_to_generate, number_of_primers_to_review, inclusion_set_threshold)
		##### NOTE: THE FLAT FILE PRODUCED BY THIS BLOCK IS def IS (self, : ):
			# <ref_seq>.scored_exclusion.InSet.<inclusion_set_threshold> 
		## APPLY PENALTY FUNCTION BY RUNNING ORACLE
		input_to_oracle = reference_sequence + ".scored_exclusion.InSet." + str(inclusion_set_threshold)
		output_from_oracle = input_to_oracle + ".oracle3penalty"
		path = utils_path
		program_call = path + 'Oracle3.py'
		os.system('python %s %s %s %s %s' %(program_call, input_to_oracle, inclusion_set_threshold, output_from_oracle, fn_blastresult))
		#### NOTE: THE FLAT FILE PRODUCED BY THIS BLOCK IS
			# <ref_seq>.scored_exclusion.InSet.<inclusion_set_threshold>.oracle3penalty
		##SORT ORACLE OUTPUT 
		#from pcr_tools import sort_by_3_priority_columns
		sorted_output_from_oracle = str(input_to_oracle + ".oracle3penalty.sorted")
		priority_column = 22 # 1 indexed : oracle penalty ##### !!! ##### HUGE CHANGE TO SIMILARITY INDEPENDENT
		second_priority_column = 12 # Inclusivity allowing for 2 mismatches
		third_priority_column = 10 # Inclusiivty wtih no mismatches 
		fourth_priority_column = 2 # 1 indexed : thermo penalty column
		In = output_from_oracle
		Out = sorted_output_from_oracle
		os.system("sort -k %i,%in -k %i,%inr -k %i,%inr -k %i,%in %s > tempsort.txt" %(priority_column, priority_column, second_priority_column, second_priority_column,third_priority_column,third_priority_column,fourth_priority_column,fourth_priority_column, In)) # SORT COMMAND, SORT FIRST ON THE priority column  numerically, next sort on the 2nd priority column 
		
		os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
		#!# ADDITIONAL OUTPUT TO SPECIAL FOLDER
		os.system('head -%i tempsort.txt > %s%s' %(number_of_primers_to_review, final_output_path, Out))
		
		#sort_by_3_priority_columns(In,Out,priority_column, second_priority_column, third_priority_column, number_of_primers_to_review, direction)
		
		## REV_SORT, FOR AN EXTENDED SET> REALLY SUPPOSED TO CONFIRM THE LESS FREQUENT. 
		revsorted_output_from_oracle = str(input_to_oracle + ".oracle3penalty.revsorted")
		In = output_from_oracle
		Out = revsorted_output_from_oracle
		#sort_by_priority_columns(output_from_oracle,revsorted_output_from_oracle,priority_column, second_priority_column, number_of_primers_to_review, 'R')
		os.system("sort -k %i,%inr -k %i,%inr -k %i,%inr -k %i,%in %s > tempsort.txt" %(priority_column, priority_column, second_priority_column, second_priority_column,third_priority_column,third_priority_column,fourth_priority_column,fourth_priority_column, In)) # SORT COMMAND, SORT FIRST ON THE priority column  numerically, next sort on the 2nd priority column 
		
		os.system('head -%i tempsort.txt > %s' %(number_of_primers_to_review, Out))
		#!#  ADDITIONAL OUTPUT TO SPECIAL FOLDER
		os.system('head -%i tempsort.txt > %s%s' %(number_of_primers_to_review, final_output_path, Out))
		
		##SELECT NON-REDUNDANT PRIMERS.
		path = utils_path
		program_call = path + 'nL_nr_assay_check.py'
		nr_in = sorted_output_from_oracle
		nr_out_specific= sorted_output_from_oracle + ".nr"
		os.system('python %s %s > %s' %(program_call, nr_in, nr_out_specific))
		nr_in = revsorted_output_from_oracle # revsorted
		nr_out_extended =  revsorted_output_from_oracle + ".nr"
		
		os.system('python %s %s > %s' %(program_call, nr_in, nr_out_extended))
		#!# ADDITIONAL OUTPUT TO SPECIAL FOLDER
		os.system('python %s %s > %s%s' %(program_call, nr_in,final_output_path, nr_out_extended))
		
		
		## FINAL OUTPUT... 
		singular_output_name = "printout.temp"
		
		
		path = utils_path
		final_output_path = str(final_output_path)
		program_call = path + 'nL_print_out.py'
		desired_primers_per_ref_seq = 8
		print_out_input = str()
		primer_type = 'specific'
		os.system('python %s %s %s > %s' %(program_call, nr_out_specific, primer_type,singular_output_name))
		os.system('head -%i %s >> %s%s'%(desired_primers_per_ref_seq, singular_output_name,final_output_path, appended_final_output_name))
		
		primer_type = 'extended'
		os.system('python %s %s %s > %s' %(program_call, nr_out_extended, primer_type,singular_output_name))
		os.system('head -%i %s >> %s%s'%(desired_primers_per_ref_seq, singular_output_name, final_output_path,appended_final_output_name))
		move_all_files_starting_with_x_to_destination(final_output_path, reference_sequence)
	fh.close()
	
	


def nL_PRIME_B(fn_nucleotide_fasta, fn_blastresult, fn_subnetwork_file, number_of_primers_to_generate, number_of_primers_to_review, utils_path, temp_path, final_output_path,  grade):
	'''
	Name: nL_PRIME_B
	
	Description: This is the primer design pipeline for those B-type reference sequences that have >= 3 inclusion set sequences but no exclusion set. We trick the system 
	
	Inputs:
		-- fn_nucleotide_fasta -- 
		-- fn_blastresult --
		-- fn_subnetwork_file --
		-- final_output_path --
		-- number_of_primers_to_generate -- 
		-- number_of_primers_to_review -- 
		-- final_output_path -- 
		-- utils_path --
		-- temp_path --  
		-- grade -- 
		
	Results :
		-- -- Major results is the primers to order file generated in the run_name labeled folder, which contains specific and extended primers that performed best and worst against the penalty function. 
	'''
	import sys
	import os 
	import datetime
	now = datetime.datetime.now()
	appended_final_output_name = 'nL_PRIMER_%s_output_for_orders_'%(grade)+ now.strftime("%Y-%m-%d_%H:%M")
	fh = open(fn_subnetwork_file, 'r')
	for line in fh:
		line = line.strip()
		units = line.split('\t')
		reference_sequence = units[0]
		print "#### REFERENCE : " +  str(reference_sequence)
		specific_exclusion_set = [] # Empty 
		specific_inclusion_set = units[6].split("|") 
		inclusion_set_threshold = units[-4]
		## DETERMINATION OF THE CANDIDATE ASSAYS THAT UNIVERSALLY MATCH THE INCLUSION SET
		from nL_Perseus3 import Perseus_part_1
		Perseus_part_1(fn_nucleotide_fasta, reference_sequence, specific_inclusion_set, number_of_primers_to_generate, number_of_primers_to_review, inclusion_set_threshold)
		##### NOTE: THE FLAT FILE PRODUCED BY THIS BLOCK IS def IS (self, : ):
			# <ref_seq>.BAE84628.1.scoredprimers.InSet.<inclusion_set_threshold> 
			# <ref_seq>.BAE84628.1.BAE84628.1.best100primers.InSet.<inclusion_set_threshold> 
		##SELECT NON-REDUNDANT PRIMERS.
		
		path = utils_path
		program_call = path + 'nL_nr_assay_check.py'
		nr_in = reference_sequence + '.best100primers.InSet.' + inclusion_set_threshold 
		nr_out_specific= nr_in + ".nr"
		os.system('python %s %s > %s' %(program_call, nr_in, nr_out_specific))
				
		## FINAL OUTPUT... 
		singular_output_name = "printout.temp"
		path = utils_path
		final_output_path = str(final_output_path)
		program_call = path + 'nL_print_out.py'
		desired_primers_per_ref_seq = 8
		print_out_input = str()
		primer_type = 'specific'
		os.system('python %s %s %s > %s' %(program_call, nr_out_specific, primer_type,singular_output_name))
		os.system('head -%i %s >> %s%s'%(desired_primers_per_ref_seq, singular_output_name,final_output_path, appended_final_output_name))
	
	fh.close()





def produce_glocks_method(fn_fasta, fn_assays_list):
	''' USED TO EFFICIENTLY SELECT GBLOCKS, GIVEN A FASTA FILE OF FULL LENGTH SEQUENCE AND A PRIMER FILE OF THE EXPECTED FORMAT LIKE THIS
	>ACF24861.1_1091_spec_F_1344	TGGCAGGCGGATAAATTCTT	>ACF24861.1_1091_spec_R_1437	CGGCACTGTCAAACCCATAA
	Returns:
		A gblock fasta file
	'''
	import sys
	from Bio import SeqIO
	handle = open(fn_fasta, "rU")
	record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()
	
	# MAKE
	fh = open(fn_assays_list, 'r')
	L1 = []
	GL = []
	Ref_L = list()
	D = dict()
	
	
	import re
	for line in fh:
		line = line.strip().split('\t')
		m0 = re.match('>([a-zA-Z0-9]+\.[0-9]+)', line[0])
		ref_seq_temp = m0.group(1)
		if ref_seq_temp not in D.keys():
			D[ref_seq_temp] = list()
			Ref_L.append(ref_seq_temp)
		#print line[0]
		m1 = re.match('>[a-zA-Z0-9]+\.[0-9]+_.*F_([0-9]+)', line[0]) #>AAC60788.1_1_spec_F_1373	
		#print line[2]
		m2 = re.match('>[a-zA-Z0-9]+\.[0-9]+_.*R_([0-9]+)', line[2]) #>AAC60788.1_1_spec_R_1457
		#print str(m1.group(1)) + "\t" + str(m2.group(1))
		D[ref_seq_temp].append([int(m1.group(1)), int(m2.group(1))])
		
		
	oh = open(fn_assays_list + ".gblocks", 'w')
	
	
	for r in Ref_L:
		print "#######%s#######" %(r)
		L1 = list(D[r])	
		print L1
		print
		captured_list = list()
		best_count = 0
		store_best = [0,0]
		for i in range(1800):
			gblock_start = i
			gblock_end = i + 490
			captured_count = 0
			for k in range(len(L1)):
				if int(L1[k][0]) - 20 > gblock_start and int(L1[k][1]) + 20 < gblock_end:
					captured_count = captured_count + 1 
			if captured_count > best_count:
				#print captured_count
				best_count = int(captured_count)
				store_best[0] = gblock_start
				store_best[1] = gblock_end
		L2 = list()
		for j in range(len(L1)):
			if int(L1[j][0]) < store_best[0] or int(L1[j][1]) > store_best[1]:
				L2.append(L1[j])
				
		print store_best
		seq = str(record_dict[r].seq)
		sub_seq = seq[store_best[0]:store_best[1]] # TAKE THE INTEVAL
		oh.write('>' + record_dict[r].description +  "_gBlock1_" + str(store_best[0]) + "_" + str(store_best[1]) + "_captures_" + str(len(L1) - len(L2)) + "\n" )
		oh.write(sub_seq + "\n")
		print str(len(L1) - len(L2)) +  " SEQUENCES WERE CAPTURED"
		
		if len(L2) > 0: # YOU FAILED TO CAPTURE ALL THE SEQUENCES IN ONE GBLOCK 
			captured_list = list()
			best_count = 0
			store_best = [0,0]
			for i in range(1800):
				gblock_start = i
				gblock_end = i + 490
				captured_count = 0
				for k in range(len(L2)):
					if int(L2[k][0]) - 20 > gblock_start and int(L2[k][1]) + 20 < gblock_end:
						captured_count = captured_count + 1 
				if captured_count > best_count:
					#print captured_count
					best_count = int(captured_count)
					store_best[0] = gblock_start
					store_best[1] = gblock_end
					
			L3 = list()
			for j in range(len(L2)):
				if int(L2[j][0]) < store_best[0] or int(L2[j][1]) > store_best[1]:
					L3.append(L2[j])
					
			print store_best
			print str(len(L2) - len(L3)) +  " SEQUENCES WERE CAPTURED"
			seq = str(record_dict[r].seq)
			sub_seq = seq[store_best[0]:store_best[1]] # TAKE THE INTEVAL
			oh.write('>' + record_dict[r].description +  "_gBlock2_" + str(store_best[0]) + "_" + str(store_best[1]) + "_captures_" + str(len(L2) - len(L3)) + "\n" )
			oh.write(sub_seq + "\n")
			
			
			
			if len(L3) > 0: # YOU FAILED TO CAPTURE ALL THE SEQUENCES IN ONE GBLOCK 
				captured_list = list()
				best_count = 0
				store_best = [0,0]
				for i in range(1800):
					gblock_start = i
					gblock_end = i + 490
					captured_count = 0
					for k in range(len(L3)):
						if int(L3[k][0]) - 20 > gblock_start and int(L3[k][1]) + 20 < gblock_end:
							captured_count = captured_count + 1 
					if captured_count > best_count:
						#print captured_count
						best_count = int(captured_count)
						store_best[0] = gblock_start
						store_best[1] = gblock_end
					
				non_captured_list = list()
				for j in range(len(L3)):
					if int(L3[j][0]) < store_best[0] or int(L3[j][1]) > store_best[1]:
						non_captured_list.append(L3[j])
			
				print store_best
				print str(len(L3) - len(non_captured_list)) +  " SEQUENCES WERE CAPTURED"
				seq = str(record_dict[r].seq)
				sub_seq = seq[store_best[0]:store_best[1]] # TAKE THE INTEVAL
				oh.write('>' + record_dict[r].description +  "_gBlock3_" + str(store_best[0]) + "_" + str(store_best[1]) + "_captures_" + str(len(L3) - len(non_captured_list)) + "\n" )
				oh.write(sub_seq + "\n")
			


def prep_gblock_IDT_orders(fn_in,start):
	import re
	import sys
	fh = open(fn_in, 'r')
	oh = open(fn_in + ".idt_ready_gblocks", "w")
	count = start
	while True:
		line = fh.readline().strip()
		if line.startswith(">"):
			count = count + 1
			m = re.match('.*_gBlock[0-9]_([0-9]+)_([0-9]+)_captures_([0-9]+)', line)
			units = line.split(" ")
			compact_name = "%s_%s:%s_C(%s)" %(units[0], m.group(1), m.group(2), m.group(3))
			sys.stdout.write("GB" + str(count) + compact_name + "\t")
			oh.write("GB" + str(count) + compact_name + "\t")
			line = fh.readline().strip()
			sys.stdout.write(line + "\t" + " ".join(map(str,units[1:]))+ "\n")
			oh.write(line + "\t" + " ".join(map(str,units[1:]))+ "\n")
		elif not line:
			break
	fh.close()
	oh.close()


def prep_assays_IDT_order(well_labels, fn_in):
	fh = open(fn_in, 'r')
	oh = open(fn_in + ".idt_ready_assays", 'w')
	forward_list = []
	reverse_list = []
	for line in fh:
		line = line.strip().split()
		well = well_labels.pop(0)
		forward_list.append( (well+"_forward" , line[0], line[1]) )
		reverse_list.append( (well+"_reverse" , line[2], line[3]) )
	fh.close()
	
	for fl in forward_list:
		oh.write("\t".join(map(str,fl))+"\n")
		print "\t".join(map(str,fl))
	
	for rl in reverse_list:
		oh.write("\t".join(map(str,rl))+"\n")
		print "\t".join(map(str,rl))
	oh.close()
	


def print_HS_nL_prime_logo():
	print '	HH    HH  SS S S SS               LL          PP P P PP  RR R R RR  II  MM      MM  EE E E E' 
	print '	HH    HH  SS                      LL          PP     PP  RR     RR  II  MM M  M MM  EE'
	print '	HH H HHH  SS S S SS    NN     NN  LL          PP P P PP  RR R R RR  II  MM  M   MM  EE E E'
	print '	HH    HH         SS    NN N   NN  LL          PP         RR RR      II  MM      MM  EE'
	print '	HH    HH         SS    NN   N NN  LL          PP         RR   RR    II  MM      MM  EE'
	print '	HH    HH  SS S S SS    NN     NN  LL L L L    PP         RR     RR  II  MM      MM  EE E E E'


def print_nL_prime_logo():
	print
	print"             LL          PP P P PP  RR R R RR  II  MM      MM  EE E E E "
	print"             LL          PP     PP  RR     RR  II  MM M  M MM  EE "
	print"  NN     NN  LL          PP P P PP  RR R R RR  II  MM  M   MM  EE E E "
	print"  NN N   NN  LL          PP         RR RR      II  MM      MM  EE "
	print"  NN   N NN  LL          PP         RR   RR    II  MM      MM  EE "
	print"  NN     NN  LL L L L    PP         RR     RR  II  MM      MM  EE E E E "  
	print 
	
def print_IDT_ordering_tool():	
	print
	print"                I   DDDD    TTTTT"
	print"                I   D   D     T"
	print"                I   D    D    T"
	print"M   M   Y Y     I   D    D    T     TTTTT   OOOOO   OOOOO   L"
	print"M M M    Y      I   D   D     T       T     O   O   O   O   L"
	print"M   M    Y      I   DDDD      T       T     OOOOO   OOOOO   LLLLL"
	print



