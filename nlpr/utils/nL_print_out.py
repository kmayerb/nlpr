# THE PURPOSE IS TO TAKE OUTPUT <ref_seq>.scored_exclusion.InSet.<inclusion_set_threshold>.oracle3penalty.revsorted
#<ref_seq>.scored_exclusion.InSet.<inclusion_set_threshold>.oracle3penalty.sorted

import sys

fn_in = sys.argv[1]
fh = open(fn_in, 'r')
#x = int(sys.argv[2])
primer_type = str(sys.argv[2])

for line in fh:
	line = line.strip() 
	name, penalty, p1, p2 , notes , start, stop = line.split()[0], line.split()[1],line.split()[2], line.split()[3], line.split()[11], line.split()[6], line.split()[7]
	#print line.split()
	if primer_type == 'specific':
		sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(name + "_spec" + "_F_" + start , p1, name +"_spec"+ "_R_" + stop , p2, penalty, line))
	elif primer_type == 'extended':
		sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(name + "_ext" + "_F_" + start , p1, name +"_ext"+ "_R_" + stop , p2, penalty, line))		
fh.close()




# 
# fh = open(fn_in, 'r')
# counter = 0
# for line in fh:
# 	line = line.strip()
# 	print line.split()
# 	name, p1, p2, notes , start, stop = line.split()[0], line.split()[2], line.split()[3], line.split()[11], line.split()[6], line.split()[7]
# #print p1
# #print p2
# 	print "%s\t%s\t%s" %(name +"_pr"+ "_R_" + stop , p2, line)
# fh.close()
# 
# 

# counter = 0
# while counter < x:
# 	line = fh.readline()
# 	counter = counter + 1
# 	line = line.strip() 
# 	name, p1, p2 , notes , start, stop = line.split()[0], line.split()[2], line.split()[3], line.split()[11], line.split()[6], line.split()[7]
# 	
# 	print "%s\t%s\t%s" %(name + "_F_" + start , p1, line)
# 
# fh.close()
# fh = open(fn_in, 'r')
# counter = 0
# while counter < x:
# 	line = fh.readline()
# 	counter = counter + 1
# 	line = line.strip()
# 	name, p1, p2 , notes , start, stop = line.split()[0], line.split()[2], line.split()[3], line.split()[11], line.split()[6], line.split()[7]
# 	
# 	print "%s\t%s\t%s" %(name + "_R_" + stop , p1, line)
# fh.close()	
	


