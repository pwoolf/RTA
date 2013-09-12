#simple script to generate synthetic data
# subject, outcome, task1_value, task2_value...
# each row is the next time step
#oucome can be 1 or 0 (redo or not)
# task_values can take any finite number of states
# 
import random , sys

def generate_random_task_series(total_task_time, task_arities, subject_id, outcome):
	task_series=[]
	for event in xrange(total_task_time):
		temp=[]
		temp.append(subject_id)
		temp.append(outcome)
		for i in task_arities:
			temp.append(random.randint(1,i))
		task_series.append(temp)
	return task_series

def impose_pattern1(outcome):
	#variable 2=1, variable 5=2 or 3
	#variable 2=2, variable 5=2 
	#gap of 1 or 2
	#variable 2=1, varaible 5=1
	max_pattern_length=5
	start=random.randint(0,len(outcome)-max_pattern_length)
	for i in xrange(max_pattern_length):
		if i==0:
			outcome[start+i][2]=1
			if random.random()>0.5:
				outcome[start+i][5]=2
			else:
				outcome[start+i][5]=3
		elif i==1:
			outcome[start+i][2]=2
			outcome[start+i][5]=2
		elif i==4:
			if random.random()>0.5:
				outcome[start+i][2]=1
				outcome[start+i][5]=1
			else:
				outcome[start+i-1][2]=1
				outcome[start+i-1][5]=1
	return outcome 

def translate_outcome_to_string(this_outcome):
	temp=[]
	for o in this_outcome:
		temp.append(",".join(map(str,o)))
	return "\n".join(temp)

def main(argv):
	task_arities		=[2]*3+[3]*1 #number of states these can take on
	num_subjects		=30 #in each group
	total_task_time		=20 #number of steps in the total task
	pattern_probability	=1.0 #probability that pattern will be introduced into a positive task series
	outcomes=[]
	#negative group
	outcome=0
	for subject_id in xrange(num_subjects):
		this_outcome=generate_random_task_series(total_task_time, task_arities, subject_id, outcome)
		outcomes.append(translate_outcome_to_string(this_outcome))
		
	#positive group:
	outcome=1
	for subject_id in xrange(num_subjects,num_subjects*2):
		this_outcome=generate_random_task_series(total_task_time, task_arities, subject_id, outcome)	
		if random.random()<pattern_probability:
			this_outcome=impose_pattern1(this_outcome)
		outcomes.append(translate_outcome_to_string(this_outcome))
	f=open("synthetic_data.txt","w")
	f.write("\n".join(outcomes))
	f.close()
	return 

if __name__ == "__main__":
	main(sys.argv[1:])
	
	




