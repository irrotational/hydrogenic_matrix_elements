import sympy.physics.hydrogen as hydrogenic
import numpy as np

def generate_lookup_tables(num_grid_points,n_max,charge,unique_nl_values,R_c=None):

	# Exponential grid construction with trapezoidal integration strips

	if R_c is None:
		R_c = 10 + (2.5) * n_max**2 # automatically calculated cutoff that seems to work well

	print('Using %d grid points with a cutoff of %.3f Bohr. Nuclear charge Z = %.3f' % (num_grid_points,R_c,charge))

	#grid_points_exp = [ ( 1e-6 * (radial_cutoff/1e-6)**(n/(num_grid_points)) - 1e-6 ) for n in range(0,num_grid_points) ] # Adapted exponential grid from: https://math.nist.gov/DFTdata/atomdata/node6.html
	grid_points_exp = [ ( R_c**(n/num_grid_points) - 1 ) for n in range(0,num_grid_points) ] # But this seems to work better
	grid_points_exp_spacing = [ (grid_points_exp[n+1]-grid_points_exp[n-1])/2 for n in range(1,num_grid_points-1) ] # spacing for nth grid point, excluding n=0 and n=num_grid_points-1 (done on next line)
	grid_points_exp_spacing.insert(0,(grid_points_exp[1]-grid_points_exp[0])/2) # n=0 spacing
	grid_points_exp_spacing.append((grid_points_exp[num_grid_points-1]-grid_points_exp[num_grid_points-2])/2) # n=num_grid_points-1 spacing

	radial_grid = grid_points_exp

	# Construct radial lookup table for the hydrogenic R_nl

	print('Initialising radial lookup tables...')
	radial_lookup_table = np.zeros([n_max,n_max,len(radial_grid)]) # n_cutoff * n_cutoff arrays each consisting of len(radial_grid) zeros. Obviously, this table is `triangular' (l_max<n).
	for n in range(1,n_max+1):
		for l in range(n):
			if (n,l) in unique_nl_values:
				for r_i, r in enumerate(radial_grid):
					radial_lookup_table[n-1][l][r_i] = hydrogenic.R_nl(n,l,r,charge)
	print('Finished initialising radial lookup tables.')

	return(radial_grid,radial_lookup_table,grid_points_exp_spacing)