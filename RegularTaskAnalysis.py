import sys, time 
try:
	import re2 as re
except ImportError:
	try:
		print "re2 not loaded, trying to install.."
		#assumes EC2 architecture, umbutu.  seems not to work on the ec2 user.
		os.system("sudo apt-get update -y")
		os.system("sudo apt-get install mercurial -y")
		os.system("sudo apt-get install python-dev -y ")
		os.system("sudo apt-get install build-essential -y ")
		os.system("sudo apt-get install g++ -y")
		os.system("sudo apt-get install git -y")
		os.system("sudo apt-get install make -y")
		os.system("sudo apt-get install gcc -y")
	
		os.system("hg clone https://re2.googlecode.com/hg re2")
		os.system('sudo make --directory="re2" test')
		os.system('sudo make --directory="re2" install')
		os.system('sudo make --directory="re2" testinstall')
		os.system("git clone git://github.com/axiak/pyre2.git")
		current_directory=os.getcwd()
		os.chdir(current_directory+"/pyre2")
		os.system("sudo python setup.py install")
		os.chdir(current_directory)
		import re2 as re
	except:
		print "Can't load re2.. using slower re."
		import re
	
		print "re2 not available, falling back to re.  "
		print "You might consider installing re2 and pyre2"
		print "\thttps://code.google.com/p/re2/"
		print "\thttps://github.com/facebook/pyre2"
 
import itertools
from multiprocessing import Pool 
import math
import getopt
from hypergeometric import log_getFETprob

def getLogPvalue(p, P, n, N):
	""" Returns a Fisher's exact pvalue 
	
	:param p: Number of positive successes
	:type p: integer
	:param P: Total number of positive cases
	:type P: integer
	:param n: Number of negative successes
	:type n: integer
	:param N: Total number of negative cases
	:type N: integer
	
	:returns: P-value  
	:rtype: float
	"""
	if (p * N > n * P):
		# apply Fisher Exact test (hypergeometric p-value)
		return math.exp(log_getFETprob(N-n, n, P-p, p)[4]);
		
	else:
		return 1.0          # pvalue = 1

def get_Bonferroni_correction_threshold(num_tests,alpha=0.05):
	"""Calculates the Bonferroni p-value correction threshold for a given number of tests at a given confidence, alpha.  
	Note that this generates a conservative false discovery estimate.
	
	:param num_tests: Total number of regular expressions tested
	:type num_tests: integer
	:param alpha: False discovery rate for threshold, default value of 0.05 
	:type alpha: float
	
	:returns: P-value threshold
	:rtype: float
	"""	
	return alpha/num_tests

	
def get_letter_set():
	"""Creates an optimized list of 241 letter symbols that starts with "a", "b", "c" .. and excludes reserved characters for regular expressions.
	
	:returns: filtered list of characters
	:rtype: list of strings
	"""
	base_map="""abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"' ,=_-!@#%&<>/`~:;"""
	test=[]
	for i in range(256):
		test.append(chr(i))			

	disallowed=".^$*+?{}[]\\|()"+ base_map
	for i in disallowed:
		test.remove(i)
	test=list(base_map)+test
	return test


def create_dataset_for_regex(data,  var_subset):
	"""Subsamples a variable subset from the raw data list and constructs a datafile in the form of a list of strings.
	
	:param data: A list of lists of strings of the form  subject_id, outcome, task_var1, task_var2..
	:type encoder: list of lists of strings
	:param var_subset: A tuple of one or more indicies indicating which variables to sub-select 
	:type var_subset: tuple of integers
	:returns: tuple of pos (list of strings), neg (list of strings), and encoder (dictionary)
	:rtype: tuple of list, list, and dictionary
	"""
	#create dictionary to see what combinations are possible
	encoder={}
	letter_set=get_letter_set()
	pos=[]
	neg=[]
	subject_id=data[0][0]
	oucome=data[0][1]
	temp=[]
	for d in data:
		obs_subject=d[0]
		if obs_subject!=subject_id:
			if outcome=="0": #neg
				neg.append("".join(temp))
			else:
				pos.append("".join(temp))
			temp=[]
			subject_id=obs_subject
		outcome=d[1]
		observed_state=[]
		for v in var_subset:
			observed_state.append(d[v])
		if not encoder.has_key(tuple(observed_state)):
			encoder[tuple(observed_state)]=letter_set[len(encoder.keys())]
		temp.append(encoder[tuple(observed_state)])
	#get the last one
	if temp:
		if outcome=="0": #neg
			neg.append("".join(temp))
		else:
			pos.append("".join(temp))
	#print pos	
	#print var_subset
	#print observed_states
	#now assign these to a letter
	return pos, neg, encoder 
	
def read_data(f_name="synthetic_data.txt"):
	"""Reads in comma delimited raw data of the form subject_id, outcome, task_var1, task_var2...  
	Comments indicated by a leading "#" are ignored.
	
	:param f_name: file name for the raw data
	:type f_name: string
	:returns: Split list of data
	:rtype: list of list of strings
	"""
	#returns a 2D array of 
	#
	# rows are time
	f=open(f_name).read().replace('\r\n', '\n').replace('\r', '\n').split('\n')
	temp=[]
	for l in f:
		if l[0]!="#":
			temp.append(l.split(","))
	return temp 
	
def get_test_patterns(encoder):
	"""Generates a list of regular expression patterns based on an encoder's alphabet
	
	:param encoder: dictionary of letter to task pattern associations
	:type encoder: dictionary
	:returns: list of regular expression patterns
	:rtype: list of strings
	"""
	alphabet=[]
	for item in encoder.keys():
		alphabet.append(encoder[item])
	#generate the basic test patterns
	degen_pair=[]
	
	gaps=["",".","..","..."]
	for i in xrange(1,2):
		for j in xrange(2):
			gaps.append(".{"+str(i)+","+str(i+2+j)+"}") #variable gaps of 2
			gaps.append(".{"+str(i)+","+str(i+1+j)+"}") #variable gaps of 1
	gaps=list(set(gaps))
	#single letter
	for a1 in alphabet:
		for a2 in alphabet:
			if a1>a2:
				degen_pair.append("["+a1+a2+"]")
			
	#print "alphabet length:{}".format(len(alphabet))
	#print alphabet
	alphabet=list(set(alphabet+degen_pair))
	#print len(alphabet)
	patterns=[]
	for a1 in alphabet:
		patterns.append(a1)
		for a2 in alphabet:
			patterns.append(a1+a2)
			for g in gaps:
				patterns.append(a1+g+a2)
			for a3 in alphabet:
				patterns.append(a1+a2+a3)
				for g in gaps:
					patterns.append(a1+g+a2+a3)
					patterns.append(a1+a2+g+a3)
	patterns=list(set(patterns))
	#print "Total patterns to test: {}".format(len(patterns))
	return patterns
	
def score_pattern(pattern, pos,neg):
	"""Calculates a one tailed p-value for a pattern using a Fisher's exact test
	
	:param pattern: regular expression search pattern
	:type pattern: string
	:param pos: list of strings where the outcome was positive
	:type pos: list of strings
	:param neg: list of strings where the outcome was negative
	:type neg: list of strings
	:returns: one tailed p-value
	:rtype: float
	"""
	p=re.compile(pattern)
	pp=0
	pn=0
	for l in pos:
		if p.search(l):
			pp+=1
	for l in neg:
		if p.search(l):
			pn+=1
	s1=getLogPvalue(pp, len(pos), pn, len(neg))
	s2=getLogPvalue(pn, len(neg),pp, len(pos))

	return min(s1,s2)
	
def get_top_results(combos,top_keep):
	""" Gets the top_keep scoring patterns from the dictionary combos
	
	:param combos: dictionary of regular expressions and scores
	:type combos: dictionary
	:param top_keep: number of top results to retain
	:type top_keep: integer
	:returns: a sorted list of the top patterns, their scores, entropies, and the encoder for each.
	:rtype: list of lists
	"""
	s=[]
	for c in combos.keys(): 
		s.append(combos[c]["s"])
	s=list(set(s))
	s.sort()
	out=[]
	for score in s[:top_keep]:
		for c in combos.keys():
			if score==combos[c]["s"]:
				temp=[]
				temp.append(score)
				temp.append(c[1])
				temp.append(c[0])
				
				temp.append(combos[c]["encoder"])
				out.append("\t".join(map(str,temp)))
	return "\n".join(out)

def trim_combos(combos,top_keep):
	""" Takes a list of regular expression combinations, and trims that list in place.
	
	:param combos: dictionary of regular expressions and associated data
	:type combos: dictionary
	:param top_keep: number of top scores to keep
	:type top_keep: integer
	:returns: dictionary of the same form as combos, just with fewer entries.
	:rtype: dictionary
	"""
	scores=[]
	for c in combos.keys():
		scores.append(combos[c]["s"])
	scores=list(set(scores))
	scores.sort()
	scores=scores[:top_keep]
	min_score=scores[-1]
	scores=set(scores)
	for c in combos.keys():
		if combos[c]["s"] not in scores:
			combos.pop(c)
	return combos, min_score

def run_pattern_search(q):
	""" core slave thread for running the analysis
	
	:param q:  a list of the inputs to test
	:type q: list
	:returns: a dictionary of regular expression combinations and their characterizations.
	:rtype: dictionary
	"""
	var_subset, data, s_min, top_keep=q
	combos={}
	pos, neg, encoder=create_dataset_for_regex(data,  var_subset)
	#print encoder
	test_patterns=get_test_patterns(encoder)
	print "Variable subset: {}\t Total patterns to test: {}".format(var_subset,len(test_patterns))
	counts=0
	for test_pattern in test_patterns:
		s=score_pattern(test_pattern, pos,neg)
		if s<s_min:
			name=(test_pattern, var_subset)
			combos[name]={}
			combos[name]["s"]=s
			combos[name]["encoder"]=encoder
			counts+=1
			if counts>top_keep: #keeps memory in check and increases the min score
				counts=0
				combos,s_min=trim_combos(combos,top_keep)
			#print "{}\t{}\t{}".format(s,var_subset,test_pattern)
	if len(combos.keys())>top_keep:
		combos,s_min=trim_combos(combos,top_keep)
	return combos

def get_combinations_of_task_indicies(num_tasks):
	""" Generates a list of all unique singe way and pairwise combinations of task indicies
	
	:param num_tasks: the total number of tasks
	:type num_tasks: integer 
	:returns: 
	:rtype:list
	"""
	task_indicies=[]
	#single way and pairs
	for j in reversed(xrange(1,3)):
		for i in itertools.combinations(range(2,num_tasks+2),j):
			task_indicies.append(i)
	return task_indicies

def count_number_of_hypotheses(task_indicies, data):
	""" Counts the total number of hypotheses to be tested given the dataset.  
		Used for calculating the multiple hypothesis correction, and allows what may be a long
		calculation to be stopped if it can't reach significance.
	
	:param num_tasks: the total number of tasks
	:type num_tasks: integer 
	:param data: the processed data
	:type data: array
	:returns: total number of hypotheses , and best case p-value
	:rtype: integer, float
	"""
	#not the most efficient way to do this, but okay.
	hypothesis_count=0
	for var_subset in task_indicies:
		pos, neg, encoder=create_dataset_for_regex(data,  var_subset)
		test_patterns=get_test_patterns(encoder)
		hypothesis_count+=len(test_patterns)
	s1=getLogPvalue(len(pos), len(pos), 0, len(neg))
	s2=getLogPvalue(len(neg), len(neg), 0, len(pos))
	
	best_p=min(s1,s2)
	return hypothesis_count, best_p


def main(argv):
	""" Main loop
	
	:param argv: input with data file
	:type argv:  list
	:returns: None
	"""
	#read in dataset
	try:
		opts, args = getopt.getopt(argv,"f:",["file="])
	except getopt.GetoptError:
		print "Need to specify an input file.  For example, for the synthetic data case:"
		print '\tpython RegularTaskAnalysis.py -f synthetic_data.txt'
		print "Or for the surgical case:"
		print '\tpython RegularTaskAnalysis.py -f Artery.reformatted.txt'
		sys.exit(2)
	input_file=""
	for opt, arg in opts:
	  if opt == '-f':
		input_file=arg 
	  	
	if not input_file:
		print "Need to specify an input file.  For example, for the synthetic data case:"
		print '\tpython RegularTaskAnalysis.py -f synthetic_data.txt'
		print "Or for the surgical case:"
		print '\tpython RegularTaskAnalysis.py -f Artery.reformatted.txt'
		sys.exit()
	try:
		data=read_data(input_file)
	except:
		print "*"*10+"Error in reading input file {}".format(input_file)+"*"*10
		print "Please check that the file is present and that it is properly formatted as a comma delimited  file.   For example:"
		print "\t0,0,1,1,2,3"
		print "\t0,0,2,2,2,1"
		print "The first entry is the subject or series ID, the rest describe the state of each task type."
		print "Lines starting with # are ignored."
		print "Entries can be anything, as long as the naming is consistent and does not contain commas in the entry itself."
		print "*"*40
		sys.exit()
	top_keep=10
	alpha=0.05
	#generate variable indicies to search
	num_tasks=len(data[0])-2
	
	task_indicies=get_combinations_of_task_indicies(num_tasks)
	
	num_tests, best_possible_p= count_number_of_hypotheses(task_indicies, data)
	bonf_threshold=get_Bonferroni_correction_threshold(num_tests,alpha)

	out_head="Minimum p-value based on Bonferroni correction of on {} tests with an alpha value of {}: {}\n".format(num_tests,alpha, bonf_threshold)
	print out_head[:-1]
	out_head2="Best possible p-value given these data: p-value: {}\n".format(best_possible_p)
	print out_head2[:-1]
	if best_possible_p>bonf_threshold:
		print "*"*30
		print "WARNING: Sample size too small to achieve signficance based on Bonferroni correction. Continuing with calculation to find the best results, but you have been warned!"
		print "*"*30
		out_head2+="WARNING: Sample size too small to achieve signficance based on Bonferroni correction.\n"
	#s_min=bonf_threshold
	s_min=alpha
	combos={}
	child_args=[]
	for var_subset in task_indicies:
		child_args.append([var_subset, data, s_min, top_keep])
	print "Testing {} variable combinations".format(len(child_args))
	
	pool=Pool()
	for o in pool.map(run_pattern_search,child_args):
		combos.update(o)
	
	output=out_head+out_head2+get_top_results(combos,top_keep)
	#print output
	print "saving results to file: {}.out.txt".format(input_file)
	f=open("{}.out.txt".format(input_file),"w")
	f.write(output)
	f.close()
	return 
	
if __name__ == "__main__":
	t0=time.time()
	main(sys.argv[1:])
	print "Total time: "+str(time.time()-t0)
