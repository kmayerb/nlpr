# ORACLE 3
# *kmb Sept 18, 2012, with the new modified scoring system
# *kmb Sept 20, 2012
# *kmb Sept 21, 2012

# >BAE84628.1_0	0.054721	GGCTATTATGCAGCGCCGTG	CAGAACTCGCGTACCCCGAA	63.949	64.004	1062	1258	197	1.0	BAE84628.1|CAD28792.1	1.0	BAE84628.1|CAD28792.1	1.0	BAE84628.1|CAD28792.1	1.0	BAE84628.1|CAD28792.1	1.0	BAE84628.1|CAD28792.1	0.005	AAO60101.1|CAD28790.2	0.047	CAR57926.1|ACH87597.1|CAR57937.1|BAF57046.1|CAR57933.1|CAR57929.1|CAJ75430.1|CAR57935.1|AAO60101.1|CAD28790.2|CAR57927.1|ACH87594.1|CAR57932.1|CAR57936.1|CAJ75435.1|CAR57934.1|CAR57930.1|CAR57931.1	0.052	CAR57926.1|ACH87597.1|ACH87599.1|CAR57937.1|BAF57046.1|CAR57933.1|CAR57929.1|CAJ75430.1|CAR57935.1|AAO60101.1|CAD28790.2|CAR57927.1|ACH87594.1|CAR57932.1|ACH87598.1|CAR57936.1|CAJ75435.1|CAR57934.1|CAR57930.1|CAR57931.1	0.054


import sys
import os 
fh = open(sys.argv[1], 'r')
inclusion_set_threshold = int(100*float(sys.argv[2])) # March 11, 13 100*input value 0.97 == 97 which is what Oracle 3 wants
oh = open(sys.argv[3], 'w')

for full_line in fh:
	#penalty_type = "weighted" # OUTDATED CODE LINE (?)
	# inclusivity = full_line.split()[9] # This is the inclusivity of the desired sequences
	fn_blastout = sys.argv[4]#'PF13486_Full_Length350to700AA.ncbi.aa.gbk.faa.blastresult' # HARDCODED, MUST BE PRESENT IN WORKING DIRECTORY!!
	full_line = full_line.strip()
	thermo_penalty = full_line.split()[1] # This is the primer3 penalty 
	# THIS BLOCK GRABS INFO ON THE ASSAY
	primary_seq = full_line.split()[0] # >BAE84628.1_0	
	assay_id = primary_seq.replace(">","") #BAE84628.1_0 (THE SPECIFIC ASSAY)	
	primary_seq = primary_seq.split("_")[0].replace(">","") #BAE84628.1 (THE SPECIFIC PARENT SEQUENCE)

	# THIS BLOCK GRABS EVERY OTHER ELEMENT 
	match_info = full_line.split()[9::2] # NUMERICAL INFORMATION ON PERCENTAGE OF HITS 
	match_hits = full_line.split()[10::2] # LIST TYPE INFORMATION ON THOSE SEQUENCES ACTUALLY HIT
	exclusion_hits= match_hits[5:] # FOCUSING ON THE EXCLUSION LIST
	#print "exclusion_hits" 
	#print exclusion_hits
	
	
	from pcr_tools import list_partition 
	# ARG [1]: LIST OF DELINEATED ITEMS
	# ARG [2]: DELINEATOR
	# ARG [3]: NAME FOR NOT HITS
	
	# RETURNS A LIST OF LISTS
	# ["A|B|C","A|B|C|D","A|B|C|D|E","A|B|C|D|E|F"] =>[['A', 'C', 'B'], ['D'], ['E'], ['F']]
	non_redundant_hits = list_partition(exclusion_hits, "|", "No_Hits")

	# DICTIONARY HOLDS THE SIMILARITY FOR EACH HIT RELATIVE TO THE PRIMARY SEQ
	# SO THAT THE SIMARITIES CAN BE EASILY LOOKED UP 
	# [['A', 'C', 'B'], ['D'], ['E'], ['F']]  =>  [['Apple', 'Cat', 'Bacon'], ['Dog'], ['Elephant'], ['Frog']]  or  [['99.9', '89.1', '70.5'], ['88.1'], ['90.1'], ['70.5']]  

	D = {}
	from pcr_tools import grab_sub_blastout
	grab_sub_blastout(primary_seq, fn_blastout, 'BLASTout.temp') # Grab the sub-section of the blast that you need
	fh = open('BLASTout.temp', 'r')
	for line in fh:
		line = line.strip() # AAC60788.1	AAC60788.1	100.00	501	0	0	1	501	1	501	0.0	1046
		subject = line.split()[1] #AAC60788.1
		similarity = float(line.split()[2]) #100.00
		D[subject] = similarity #{AAC60788.1: 100.00}
	fh.close()
	#print D # ERROR CATCHING July 28. 2014

	
	def list_swap(L):
	# A SIMPLE FUNCTION, GIVEN ANY LIST LOOK AND REPLACE EACH ELEMENT WITH ITS DICTIONARY VALUE. IF NO DICTIONARY VALUE EXISTS, REPLACE WITH ZERO 
		OL = []
		for l in L:
			try:
				f = D[l]
				#s = "%.2f" %(f)
				OL.append(f)
			except KeyError:
				OL.append(0) # IF YOU GET A KEY ERROR YOU ARE HITTING A SEQUENCE THAT YOU SHOULD NOT SINCE IT IS NOT EVEN IN THE BLAST RESULT FOR THE PRIMARY SEQ
		return OL

	# THIS APPLIES THIS FUNCTION TO A LIST OF LISTS 
	non_redundant_hits_similarity = [list_swap(j) for j in non_redundant_hits]
	
	
	
	#print "non_redundant_hits"
	#print non_redundant_hits 
    
	#print "non_redundant_hits_similarity "
	#print non_redundant_hits_similarity 
    	
	k = 4 
	def list_difference(L):
		return [(inclusion_set_threshold-l)**k for l in L]
	non_redundant_hits_difference = [list_difference(j) for j in non_redundant_hits_similarity]	

	# Assign penalties 
	OL = []
	penalty = [1.0, 0.9, 0.8, 0.5, 0.25]
	for x in range(5):
		L = non_redundant_hits_difference[x]
		p = penalty[x]
		OL.append([p*l for l in L])
	#print "non_redundant_hits_difference"
	#print non_redundant_hits_difference
	#print "OL"
	#print OL

	final_list = [item for sublist in OL for item in sublist] 
	final_penalty = sum(final_list)
	#final_worst_penalty = max(final_list)
	#print "final_penalty"
	#print final_penalty	# Sum penalties 
	#print int(final_penalty)
	oh.write(full_line + "\t" + str(final_penalty) + "\n")
fh.close()
oh.close()



# 
# inclusion_hit_list = full_line.split()[12].split("|") #LIST
# exclusion_hit_list = full_line.split()[-1].split("|") #LIST
# 
# D = {} # DICTIONARY HOLDS THE SIMILARITY FOR EACH
# from pcr_tools import grab_sub_blastout
# grab_sub_blastout(primary_seq, fn_blastout, 'BLASTout.temp') # Grab the sub-section of the blast that you need, next you will add a column with hits and non hits
# fh = open('BLASTout.temp', 'r')
# for line in fh:
# 	line = line.strip() # AAC60788.1	AAC60788.1	100.00	501	0	0	1	501	1	501	0.0	1046
# 	subject = line.split()[1] #AAC60788.1
# 	similarity = float(line.split()[2]) #100.00
# 	D[subject] = similarity #{AAC60788.1: 100.00}
# fh.close()
# 
# red_card = 0
# penalty_list = []
# for e in exclusion_hit_list:
# 	try:
# 		similarity_score = D[e]
# 		penalty_list.append(similarity_score)
# 	except KeyError:
# 		sys.stderr.write("For assay %s DID NOT FIND %s IN BLAST RESULT!!!!!!\n" %(assay_id, e))
# 		red_card = 1 
# 		penalty_list.append(0)
# #print penalty_list		
# 
# if red_card == 1:
# 	continue
# elif penalty_type is "simple":
# 	k = 2
# 	worst_offender = min(penalty_list)
# 	final_penalty = (100 - worst_offender)**k
# 	#print "final_penalty ('simple'): %f" %(final_penalty)
# 	oh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(full_line,"assay_id",assay_id, "final_simple_penalty", final_penalty, "thermo_penalty", thermo_penalty, "inclusivity", inclusivity))
# elif penalty_type is "cumulative":
# 	k = 2
# 	final_penalty = 0 
# 	for p in penalty_list:
# 		final_penalty += ((100 - p)**k) 
# 	#print "final_penalty ('cumulative'):: %f" %(final_penalty)
# 	oh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(full_line,"assay_id",assay_id, "final_cumulative_penalty", final_penalty, "thermo_penalty", thermo_penalty, "inclusivity", inclusivity))
# elif penalty_type is "weighted":
# 	k = 2
# 	final_penalty = 0 
# 	for p in penalty_list:
# 		final_penalty += ((100 - p)**k) 
# 	#print "final_penalty ('cumulative'):: %f" %(final_penalty)
# 	final_weighted_penalty = final_penalty / float(len(penalty_list))
# 	#print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(full_line,"assay_id",assay_id, "final_weighted_penalty", final_penalty, "thermo_penalty", thermo_penalty, "inclusivity", inclusivity)
# 	oh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(full_line,"assay_id",assay_id, "final_weighted_penalty", final_weighted_penalty, "thermo_penalty", thermo_penalty, "inclusivity", inclusivity))


	
	
	


	
