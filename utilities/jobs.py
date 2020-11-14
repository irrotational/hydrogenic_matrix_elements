###########################################################################################
# Parse the user's job list (either specified with argument -todo or a file of jobs supplied with -f)

def create_job_list(file=None,todo_list=None,n_max=None,only_these_l=None,only_these_m=None):

	job_list = [] # This will be a list of shape (num_elements_to_calculate,4,3)

	if file:
		todo_list = []
		job_list_file = open(file).readlines()
		for element in job_list_file:
			todo_list.append(element)

	if todo_list:
		for element in todo_list:
			element = element.split(" ")
			element = [ x.strip() for x in element ]
			el = []
			for x in element:
				x=x.strip("(")
				x=x.strip(")")
				x=x.split(",")
				y = []
				for z in x:
					y.append( int(z) )
				y = tuple(y)
				el.append(y)
			job_list.append(el)

	# If user did not specify with -todo or a file using -f, generate a job_list...

	if (not todo_list) and (not file):

		# If no l or m restrictions are specified, the user wants to use the full basis set up to n_max.
		if only_these_l is None:
			only_these_l = range(0,n_max)
		if only_these_m is None:
			only_these_m = range(-(n_max-1),n_max)

		basis_set_size=0
		for j in range(1,n_max+1): # calculate basis_set_size
			num = j**2 # For each n, we have 2 * n^2 orbitals
			basis_set_size += num

		state_list = [0] * basis_set_size
		count = 0 # Building the state list
		for n in range(1,n_max+1):
			for l in range(0,n):
				if l in only_these_l:
					for m in range(-l,l+1):
						if m in only_these_m:
							state_list[count] = (n,l,m) # Add the state to the state_list
							count += 1

		state_list = [ x for x in state_list if x != 0 ]

		job_list = []
		for orbital_1 in state_list:
			for orbital_2 in state_list:
				for orbital_3 in state_list:
					for orbital_4 in state_list:
						job_list.append( [orbital_1,orbital_2,orbital_3,orbital_4] )

	R_nls = [] # A list holding all the (n,l) pairs (used to determine the minimal set of radial integrals we need to compute)
	for job in job_list:
		for nlm in job:
			R_nls.append( nlm[0:2] )

	return(job_list,R_nls)

