# THE PURPOSE OF THIS SCRIPT IS TO SELECT ASSAYS FROM A	LIST OF PRIMER COMBINATIONS, SUCH THAT NO MORE THAN TWO OF THE ASSAYS TARGET EXACTLY APROXIMATELY THE SAME POSITIONS .
import sys
fn = sys.argv[1]
fh = open(fn, 'r')

def plus_or_minus(x,h):
	L = []
	for i in range(h):
		L.append(int(x-i))
		L.append(int(x+i))
	return list(set(L))

def lists_overlap3(a, b):
    return bool(set(a) & set(b))

forbidden_range_F = []
forbidden_range_R = []
forbidden_range_F2 = []
forbidden_range_R2= []
forbidden_range_F3 = []
forbidden_range_R3 = []

# Take the best hit
line = fh.readline()
line = line.strip()
print line 

for line in fh:
	forbidden_range_F3 = list(forbidden_range_F2)
	forbidden_range_R3 = list(forbidden_range_R2)
	forbidden_range_F2 = list(forbidden_range_F)
	forbidden_range_R2 = list(forbidden_range_R)
	#print "#####"
	#print forbidden_range_F2
	#print forbidden_range_F3
	#print "#####"
	line = line.strip()
	start = int(line.split()[6])
	end = int(line.split()[7])
	forbidden_range_F.append(start)
	forbidden_range_R.append(end)
	test_F = plus_or_minus(int(start),4)
 	test_R = plus_or_minus(int(end),4)

	if lists_overlap3(test_F, forbidden_range_F2) and lists_overlap3(test_R,forbidden_range_R2):
		pass
	else:
		print line
fh.close()
