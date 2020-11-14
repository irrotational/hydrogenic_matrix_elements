import random
from sympy.physics.wigner import wigner_3j
from jobs import create_job_list
from make_grids import generate_lookup_tables
import numpy as np
import time
import multiprocessing as mp
from multiprocessing import Manager
import argparse
import csv

CODE_START_TIME = time.time()

parser=argparse.ArgumentParser()
parser.add_argument('-write_to_file',type=str) # If supplied, will generate the dictionary of matrix elements and write it to a python file, which can be read in for future runs. Must be a directory string.
parser.add_argument('-N',type=int,default=500) # Number of grid points in the exponential grid for on-the-fly grid generation.
parser.add_argument('-radial_cutoff',type=float,default=None) # Cutoff of each R_nl in Bohr radii for on-the-fly grid generation.
parser.add_argument('-todo',type=str,nargs='+',default=None)
parser.add_argument('-num_cores',type=int) # If num_cores is supplied, multiprocessing is turned on and will be parallelised over num_cores processes. By default, it's absent (no MP).
parser.add_argument('-Z',default=1,type=float) # Nuclear charge Z. This need not be of type int.
parser.add_argument('-f',type=str) # File containing matrix elements to be calculated
parser.add_argument('-n_cutoff',type=int,default=None) # Basis set size.
parser.add_argument('-restrict_l',nargs='+',type=int) # Only basis states with l values in this space-separated list will be used.
parser.add_argument('-restrict_m',nargs='+',type=int) # Only basis states with m values in this space-separated list will be used.
args=parser.parse_args()
write_to_file=args.write_to_file
N=args.N
radial_cutoff=args.radial_cutoff
num_cores=args.num_cores
Z=args.Z
todo=args.todo
f=args.f
restrict_l=args.restrict_l
restrict_m=args.restrict_m
n_cutoff=args.n_cutoff

# Get a list of jobs (matrix elements) to compute
# job_list is a list of shape (num_elements_to_calculate,4,3)

if f: # If user supplies a job file
	job_list,R_nls = create_job_list(file=f)
elif todo: # If user supplies the matrix elements at the command line
	job_list,R_nls = create_job_list(todo_list=todo)
elif n_cutoff: # If user specifies only n_cutoff
	job_list,R_nls = create_job_list(n_max=n_cutoff,only_these_l=restrict_l,only_these_m=restrict_m)
else:
	print('Error: No matrix elements requested! Aborting.')
	exit()

n_cutoff = max( np.array(R_nls)[:,0] ) # n_cutoff is the highest n value we need to consider

print('Number of matrix elements todo: %d' % (len(job_list)))

# Now construct the radial grids and radial lookup tables
# radial_lookup_table entries are of the form: [n-1][l][r_i], so e.g. the jth entry for (n=2,l=1) would be [1][1][j]

radial_grid,radial_lookup_table,grid_spacing = generate_lookup_tables(num_grid_points=N,n_max=n_cutoff,charge=Z,unique_nl_values=R_nls,R_c=radial_cutoff)

###########################################################################################
# Functions that check for symmetries in the integrals

def get_radial_integral(n1,l1,n2,l2,n3,l3,n4,l4,l): # Retrieves a radial integral from the four_radial_dictionary above, accounting for aliases. Will return None if the integral or equivalent alias does not exist.

	perm1 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n1,l1,n2,l2,n3,l3,n4,l4,l) # Identity
	if (perm1 in four_radial_dictionary):
		return four_radial_dictionary[perm1]
	perm2 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n2,l2,n1,l1,n3,l3,n4,l4,l)
	if (perm2 in four_radial_dictionary):
		return four_radial_dictionary[perm2]
	perm3 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n1,l1,n2,l2,n4,l4,n3,l3,l)
	if (perm3 in four_radial_dictionary):
		return four_radial_dictionary[perm3]
	perm4 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n2,l2,n1,l1,n4,l4,n3,l3,l)
	if (perm4 in four_radial_dictionary):
		return four_radial_dictionary[perm4]
	perm5 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n3,l3,n4,l4,n1,l1,n2,l2,l)
	if (perm5 in four_radial_dictionary):
		return four_radial_dictionary[perm5]
	perm6 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n4,l4,n3,l3,n1,l1,n2,l2,l)
	if (perm6 in four_radial_dictionary):
		return four_radial_dictionary[perm6]
	perm7 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n3,l3,n4,l4,n2,l2,n1,l1,l)
	if (perm7 in four_radial_dictionary):
		return four_radial_dictionary[perm7]
	perm8 = '%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n4,l4,n3,l3,n2,l2,n1,l1,l)
	if (perm8 in four_radial_dictionary):
		return four_radial_dictionary[perm8]

def get_four_orbitals_matrix_element(orbital_1, orbital_2, orbital_3, orbital_4): # Retrieves a four_orbital matrix element from four_orbitals_dictionary, accounting for aliases. Again, will return None if the integral or equivalent alias does not exist.
	
	n1,n2,n3,n4 = orbital_1[0],orbital_2[0],orbital_3[0],orbital_4[0]
	l1,l2,l3,l4 = orbital_1[1],orbital_2[1],orbital_3[1],orbital_4[1]
	m1,m2,m3,m4 = orbital_1[2],orbital_2[2],orbital_3[2],orbital_4[2]

	perm1 = '%s,%s,%s,%s' % (orbital_1, orbital_2, orbital_3, orbital_4) # Identity
	if (perm1 in four_orbitals_dictionary):
		return four_orbitals_dictionary[perm1]
	perm2 = '%s,%s,%s,%s' % (orbital_2, orbital_1, orbital_4, orbital_3) # Complex conjugate
	if (perm2 in four_orbitals_dictionary):
		return four_orbitals_dictionary[perm2].conjugate()

# This function is the one that actually computes the requested matrix element

def four_orbitals_matrix_element(orbital_1, orbital_2, orbital_3, orbital_4): # Calculates < 1(r)3(r') |  1/|r-r'|  | 2(r)4(r') >

	x = get_four_orbitals_matrix_element(orbital_1, orbital_2, orbital_3, orbital_4) # First check if it has already been done

	if x is None: # If not already done, we'll have to compute it

		n1,n2,n3,n4 = orbital_1[0],orbital_2[0],orbital_3[0],orbital_4[0]
		l1,l2,l3,l4 = orbital_1[1],orbital_2[1],orbital_3[1],orbital_4[1]
		m1,m2,m3,m4 = orbital_1[2],orbital_2[2],orbital_3[2],orbital_4[2]

		total = 0

		l_range_l1_l2 = np.arange( abs(l1-l2),abs(l1+l2)+1 ) # Triangle condition l range due to l1 and l2
		l_range_l3_l4 = np.arange( abs(l3-l4),abs(l3+l4)+1 ) # Ditto for l3 and l4
		l_range = np.intersect1d(l_range_l1_l2,l_range_l3_l4) # Only include terms in the sum that are in both ranges above; others are zero
		l_range = filter(lambda x: (x+l1+l2)%2==0, l_range) # Now remove terms that don't satisfy (l + l1 + l2) is even; gives zero 3-j symbol
		l_range = filter(lambda x: (x+l3+l4)%2==0, l_range) # Ditto for terms that don't satisfy (l + l3 + l4) is even

		for l in l_range:

			radial_total = 0
			rad = get_radial_integral(n1,l1,n2,l2,n3,l3,n4,l4,l)

			if (rad is None): # If it's not in the dictionary, add it.

				for r_i, r in enumerate(radial_grid):
					if r <= 1e-100: # Integrand contribution is very small below here, and leads to division by zero => Ignore
						continue
					for r_j, rp in enumerate(radial_grid):
						if rp <= 1e-100: # Ditto
							continue

						if r >= rp:
							prefactor = ( rp**l/r**(l+1) )
						else:
							prefactor = ( r**l/rp**(l+1) )

						R_nl_part  = radial_lookup_table[n1-1][l1][r_i] # Remember the n indexing is shifted by 1 (see above)
						R_nl_part *= radial_lookup_table[n2-1][l2][r_i]
						R_nl_part *= radial_lookup_table[n3-1][l3][r_j]
						R_nl_part *= radial_lookup_table[n4-1][l4][r_j]

						radial_total += r**2 * rp**2 * grid_spacing[r_i] * grid_spacing[r_j] * R_nl_part * prefactor # for non-linear grid

				four_radial_dictionary['%d,%d,%d,%d,%d,%d,%d,%d,l=%d' % (n1,l1,n2,l2,n3,l3,n4,l4,l)] = radial_total
			else:
				radial_total = rad

			for m in range(-l,l+1):
				if (m2 - m - m1 == 0) and (m - m3 + m4 == 0): # The 3-j symbols are zero unless these are satisfied; saves some time.

					angular_part = (-1)**(m+m1+m3)
					angular_part *= np.sqrt( (2*l1+1)*(2*l2+1)*(2*l3+1)*(2*l4+1) )
					angular_part *= float( wigner_3j(l,l1,l2,0,0,0) ) # Need to convert these to float else they're a 'sympy.core.numbers' type
					angular_part *= float( wigner_3j(l,l3,l4,0,0,0) ) # Ditto
					angular_part *= float( wigner_3j(l,l1,l2,-m,-m1,m2) ) # Ditto
					angular_part *= float( wigner_3j(l,l3,l4,m,-m3,m4) ) # Ditto

					total += radial_total * angular_part

		x = total

	four_orbitals_dictionary[ '(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d)' % (n1,l1,m1,n2,l2,m2,n3,l3,m3,n4,l4,m4) ] = x
	return Z * x

# Explain notation to the user
print('\n')
print('Matrix Element Notation: (n_a,l_a,m_a),(n_b,l_b,m_b),(n_c,l_c,m_c),(n_d,l_d,m_d)')
print('\n')
print('Is Equivalent To: < a(r)c(r\') |  1/|r-r\'|  | b(r)d(r\') >')
print('\n')

###########################################################################################
# Functions for multiprocessing. These are only called if num_cores is not None.

def chunk_processes(job_list):
	random.shuffle(job_list) # Shuffles the job list (inplace) so that no one core is biased towards easier/trickier calculations
	chunk_size = int( len(job_list) / num_cores )
	chunks = [job_list[j*chunk_size:(j+1)*chunk_size+1] for j in range(num_cores-1)] # Split job_list into num_cores chunks
	chunks.append(job_list[(num_cores-1)*chunk_size:]) # The last chunk will in general have an unequal size
	return chunks

def calculate_chunk(chunk):
	for job in chunk:
		print( 'Calculating Matrix Element:', *job, '\n', 'Result:', four_orbitals_matrix_element(*job) )
		print('\n')

###########################################################################################
# Calculation of the elements

orbitals_dict_starttime = time.time()

if num_cores is None:
	four_radial_dictionary = {} # Will be populated with four-radial integrals
	four_orbitals_dictionary = {} # Will be populated with four-orbital matrix elements
	for job in job_list:
		print( 'Calculating Matrix Element:', *job )
		print( ' Result:',four_orbitals_matrix_element(*job) )
		print('\n')

else:
	print('Running calculation in parallel with %d processes\n' % (num_cores))
	manager=Manager()
	four_radial_dictionary = manager.dict() # In the multiprocessing case, we need a shared-memory dictionary type
	four_orbitals_dictionary = manager.dict() # Ditto
	matrix_element_joblist = chunk_processes(job_list)
	processes=[]
	for chunk in matrix_element_joblist:
		p = mp.Process(target=calculate_chunk, args=(chunk,))
		processes.append(p)
		p.start()

	for process in processes:
		process.join()

###########################################################################################
# Write to file, if requested

if write_to_file:
	csv_outdict = csv.writer(open('./'+write_to_file+'.csv',"w"))
	for key, val in four_orbitals_dictionary.items():
		csv_outdict.writerow([key,val])
	print('Matrix Elements have been written to %s.csv' % (write_to_file))

CODE_STOP_TIME = time.time()
print('Total Run Time: %f seconds' % (CODE_STOP_TIME-CODE_START_TIME))

