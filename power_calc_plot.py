import sys, time 
import itertools
from multiprocessing import Pool 
import math
import getopt

from hypergeometric import log_getFETprob

import matplotlib.pyplot as plt
from pylab import *

def getLogPvalue(p, P, n, N):
	""" Return log of hypergeometric pvalue of #pos >= p
			p = positive successes
			P = positives
			n = negative successes
			N = negatives
	"""
	if (p * N > n * P):
		# apply Fisher Exact test (hypergeometric p-value)
		return math.exp(log_getFETprob(N-n, n, P-p, p)[4]);
		
	else:
		return 1.0          # pvalue = 1

def f(x):  #factorial
	return math.factorial(x)
	
		
def fish(a,b,c,d): 
	temp= math.pow(math.e,math.log(f(a+b))+math.log(f(a+c))+math.log(f(c+d))+math.log(f(b+d))-math.log(f(a))-math.log(f(b))-math.log(f(c))-math.log(f(d))-math.log(f(a+b+c+d)))
	return temp

def fullfish(a,b,c,d): # calculates the one tailed fishers exact p-value
	# TABLE:			motif
	#				____pres____abs__
	#	  in	pos	|	a		b
	#			neg	|	c		d
	ptemp=fish(a,b,c,d)
	atemp=a+1
	btemp=b-1
	ctemp=c-1
	dtemp=d+1
	
	while btemp>=0 and ctemp>=0:
		#print "********"
		#print atemp,btemp,ctemp,dtemp
		ptemp+=fish(atemp,btemp,ctemp,dtemp) #I'm assuming more extreme here is that more varables are in the set.
		atemp+=1
		btemp-=1
		ctemp-=1
		dtemp+=1
	return ptemp	

def get_Bonferroni_correction_threshold(num_tests,alpha=0.05):
	#gets the minimum acceptance p-value threshold for a given alpha.
	# note that the Bonferroni correction is conservative
	
	return alpha/num_tests

	
def get_letter_set():
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
	#create dictionary to see what combinations are possible
	observed_states={}
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
		if not observed_states.has_key(tuple(observed_state)):
			observed_states[tuple(observed_state)]=letter_set[len(observed_states.keys())]
		temp.append(observed_states[tuple(observed_state)])
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
	return pos, neg, observed_states 
	
def read_data(f_name="synthetic_data.txt"):
	#returns a 2D array of 
	# subject_id, outcome, task_var1, task_var2...
	# rows are time
	f=open(f_name).read().replace('\r\n', '\n').replace('\r', '\n').split('\n')
	temp=[]
	for l in f:
		if l[0]!="#":
			temp.append(l.split(","))
	return temp 
	
def get_test_patterns(encoder):
	alphabet=[]
	for item in encoder.keys():
		alphabet.append(encoder[item])
	#generate a bunch of basic test patterns
	degen_pair=[]
	degen_triple=[]
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
			for a3 in alphabet:
				if a1>a2>a3:
					degen_triple.append("["+a1+a2+a3+"]")
	degen_triple=[]
	alphabet=list(set(alphabet+degen_pair+degen_triple))
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
			
	#print "Total patterns to test: {}".format(len(patterns))
	return patterns
	
def score_pattern(pattern, pos,neg):
	p=re.compile(pattern)
	pp=0
	pn=0
	for l in pos:
		if p.search(l):
			pp+=1
	for l in neg:
		if p.search(l):
			pn+=1
	#print "{}\t{}".format(pp,pn)
	s1=getLogPvalue(pp, len(pos), pn, len(neg))
	s2=getLogPvalue(pn, len(neg),pp, len(pos))
	#s1=log_getFETprob(pp,len(pos)-pp,pn,len(neg)-pn);# fullfish(pp,len(pos)-pp,pn,len(neg)-pn) #left tail, right tail
	#s2=fullfish(len(pos)-pp,pp,len(neg)-pn,pn)
	
	return min(s1,s2)
	
def get_top_results(combos,top_keep):
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
	var_subset, data, s_min, top_keep=q
	print var_subset
	
	combos={}
	pos, neg, encoder=create_dataset_for_regex(data,  var_subset)
	#print encoder
	test_patterns=get_test_patterns(encoder)
	print "total patterns to test: {}".format(len(test_patterns))
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
	task_indicies=[]
	#single way and pairs
	for j in reversed(xrange(1,3)):
		for i in itertools.combinations(range(2,num_tasks+2),j):
			task_indicies.append(i)
	return task_indicies

def count_number_of_hypotheses(task_indicies, data):
	#not the most efficient way to do this, but okay.
	hypothesis_count=0
	for var_subset in task_indicies:
		pos, neg, encoder=create_dataset_for_regex(data,  var_subset)
		test_patterns=get_test_patterns(encoder)
		hypothesis_count+=len(test_patterns)
	best_p=get_best_case_p_value(pos,neg)
	return hypothesis_count, best_p

def get_best_case_p_value(pos,neg):
	#returns the best case p-value
	s1=getLogPvalue(len(pos), len(pos), 0, len(neg))
	s2=getLogPvalue(len(neg), len(neg), 0, len(pos))
	return min(s1,s2)

def get_iso_data_for_figure(bonf_threshold):
	np_rates=[0, 0.01, 0.05, 0.10, 0.2, 0.3, 0.4, 0.5]
	#tot_samples=[i*10 for i in range(1,21)]
	tot_samples=[i*2 for i in range(1,151)]
	pp_rates=[i/100.0 for i in range(101)] #dumb and slow, but okay.
	plot_data={}
	for np in np_rates:
		plot_data[np]=[]
		
	for num_samples in tot_samples:
		pos_len=num_samples/2
		neg_len=num_samples-pos_len
		for np_rate  in np_rates:
			best_distance=-1
			best_pp_rate=-1.0
			for pp_rate in pp_rates:
				pp=int(pos_len*pp_rate)
				pn=int(neg_len*np_rate)
				s1=getLogPvalue(pp, pos_len, pn, neg_len)
				#s2=getLogPvalue(pn, neg_len,pp, pos_len)
				#s1=fullfish(pp,pos_len-pp,pn,neg_len-pn) #exact
				#s2=fullfish(pos_len-pp,pp,neg_len-pn,pn) #exact
				s2=1.00
				#s1=1.00
				score=min(s1,s2)
				if best_distance==-1:
					best_distance=(score-bonf_threshold)**2
					best_pp_rate=pp_rate 
				elif (score-bonf_threshold)**2<best_distance:
					best_pp_rate=pp_rate 
					best_distance=(score-bonf_threshold)**2
			if best_pp_rate!=0.0 and best_pp_rate!=1.0:
			
				plot_data[np_rate].append([num_samples,best_pp_rate,best_distance])
	#print plot_data[0.00]
	return plot_data
		
def draw_iso_plot(data):
	
	out=[]
	for np_rate in data.keys():
		parts=data[np_rate]
		out.append([i[0] for i in parts])
		out.append([i[1]*100 for i in parts])
		out.append('k-')
	plot(*out)	
	#plot(x1, y1 , 'k--',x2, y2 , 'k--')
	xlabel('Total number of samples')
	ylabel('% matches to positive ')
	title('Data')
	
	savefig("iso_power.pdf")
	#show()
	return 
	#first negatives
	y_pos=1
	box_size=2
	max_width=0
	for d in data:
		if data[d]["outcome"]=="0":
			x_pos=0
			
			for obs_line in data[d]["d"]:
				rel_y_pos=y_pos
				x_pos+=box_size
				for obs in obs_line:
					rel_y_pos+=box_size
					parts.append(plt.Rectangle((x_pos,rel_y_pos),box_size,box_size,color=color_map[obs]))
			y_pos=rel_y_pos+box_size
			if x_pos>max_width:
				max_width=x_pos
				
	y_pos=rel_y_pos+box_size*10			
	for d in data:
		if data[d]["outcome"]=="1":
			x_pos=0
			
			for obs_line in data[d]["d"]:
				rel_y_pos=y_pos
				x_pos+=box_size
				for obs in obs_line:
					rel_y_pos+=box_size
					parts.append(plt.Rectangle((x_pos,rel_y_pos),box_size,box_size,color=color_map[obs]))
			y_pos=rel_y_pos+box_size
			if x_pos>max_width:
				max_width=x_pos
	fig = plt.gcf()			
	
def main(argv):
	#read in dataset
	try:
		opts, args = getopt.getopt(argv,"f:",["file="])
	except getopt.GetoptError:
		print "Need to specify an input file"
		print '\tpython power_calc_plot.py -f synthetic_data.txt'
		sys.exit(2)
	input_file=""
	for opt, arg in opts:
	  if opt == '-f':
		input_file=arg 
	  	
	if not input_file:
		print "Need to specify an input file"
		print '\tpython power_calc_plot.py -f synthetic_data.txt'
		sys.exit()
	data=read_data(input_file)
	top_keep=10
	alpha=0.05
	#generate variable indicies to search
	num_tasks=len(data[0])-2
	
	task_indicies=get_combinations_of_task_indicies(num_tasks)
	
	num_tests, best_possible_p= count_number_of_hypotheses(task_indicies, data)
	print num_tests
	return 
	bonf_threshold=get_Bonferroni_correction_threshold(num_tests,alpha)

	
	print bonf_threshold
	data=get_iso_data_for_figure(bonf_threshold)
	draw_iso_plot(data)
	


if __name__ == "__main__":
	t0=time.time()
	main(sys.argv[1:])
	print "Total time: "+str(time.time()-t0)
