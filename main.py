#################################################################################
#									Imports										#
#################################################################################
############################### Internal Libraries ##############################
import os
import random
import time
from collections import deque
from textwrap import fill as wrap
import re
import sys
############################### External Libraries ##############################
try:
	import numpy as np
except ImportError:
	print('Please install mandatory module numpy\npip install numpy')
	exit(0)

try:
	from tqdm import trange, tqdm
except:
	print('Please install optional module tqdm\n'+
		  'pip install tqdm\nTo have nice progress bar while this script runs')
	tqdm	= lambda x, *args, **kwargs: x
	trange	= lambda x, *args, **kwargs: range(x)

try:
	import matplotlib.pyplot as plt
	import matplotlib
except:
	print('Please install optional module matplotlib\n'+
		  'pip install matplotlib\nTo view the plots')
#################################################################################
#								uf20_91 selection								#
#################################################################################
# random.seed(195638)
# random.sample(list(range(1,1001)), 3)
# [385, 876, 961]
#################################################################################
#							Helper function to plot RTDs						#
#################################################################################
def plot_rtd(logs, _use_file=None, _use_legend=None, _make_file=None):
	if _use_file:
		logs = np.loadtxt(_use_file, dtype=float, delimiter=',')

	if _use_legend is True:
		# Reads the searcher's legend legend from the header of the file
		# Format: '# 0: GSAT+Tabu   1: GWalkSAT'
		# Parsed: ['GSAT+Tabu', 'GWalkSAT']
		with open(_use_file, 'r') as f:
			_use_legend = [item
						for item in f.readline().strip(' #\n').split(' ')
						if (item!='' and item[0]=='G')]

	########################## Setting up the subplots ##########################
	plt.rcParams["figure.figsize"] = (23,10)

	plt.subplot(2,4,1, ylim=(0,1), xlabel='Search steps →', ylabel='P(solve) →')
	plt.grid(True, which='both')

	plt.subplot(2,4,3, ylim=(0,1), xlabel='Runtime →', ylabel='P(solve) →')
	plt.grid(True, which='both')

	#                             semi log scales                              #
	plt.subplot(2,4,2, xscale='log', ylim=(0,1),
					   xlabel='Search steps (log scale) →', ylabel='P(solve) →')
	plt.grid(True, which='both')

	plt.subplot(2,4,4, xscale='log', ylim=(0,1),
					   xlabel='Runtime (log scale) →', ylabel='P(solve) →')
	plt.grid(True, which='both')

	#                             log - log scales                             #
	plt.subplot(2,4,5, xscale='log', yscale='log',
		   xlabel='Search steps (log scale) →', ylabel='P(solve) (log scale) →')
	plt.grid(True, which='both')

	plt.subplot(2,4,7, xscale='log', yscale='log',
				xlabel='Runtime (log scale) →', ylabel='P(solve) (log scale) →')
	plt.grid(True, which='both')

	plt.subplot(2,4,6, xscale='log', yscale='log',
	   xlabel='Search steps (log scale) →', ylabel='1 - P(solve) (log scale) →')
	plt.grid(True, which='both')

	plt.subplot(2,4,8, xscale='log', yscale='log',
			xlabel='Runtime (log scale) →', ylabel='1 - P(solve) (log scale) →')
	plt.grid(True, which='both')

	# Input file columns:
	#	0		1			2			3			4				5
	# round, searcher, sol_reached, flips_taken, comparison_runs, time_taken
	for item in np.unique(logs[:,1]): # Unique search algorithms
		steps_x = np.sort(logs[:,3][(logs[:,1]==item) & (logs[:,2]==1)])
		# logs[:,3]			# flips_taken data, filtered with following slicing
		# logs[:,1]==item	# filters the data for specific algorithm
		# logs[:,2]==1		# filters the data for only solved iterations
		total_runs = len(steps_x)
		p_solve_y = np.array([i/total_runs for i in range(total_runs)])
		#						Count(iteration)
		# to get P(solve) = ----------------------------
		#					Count(All Solved Iterations)

		if _use_legend:
			if isinstance(_use_legend, list):
				_label = _use_legend[int(item)]
			else:
				_label = _use_legend[item] # For dictionary based legend
		else:
			_label = item # Use the code itself for no legend
		################################ Plotting ###############################
		plt.subplot(2,4,1)
		plt.plot(steps_x, p_solve_y, label=_label)
		plt.subplot(2,4,2)
		plt.plot(steps_x, p_solve_y, label=_label)
		plt.subplot(2,4,5)
		plt.plot(steps_x, p_solve_y, label=_label)
		plt.subplot(2,4,6)
		p_solve_y = 1 - p_solve_y
		plt.plot(steps_x, p_solve_y, label=_label)

		steps_x = np.sort(logs[:,5][(logs[:,1]==item) & (logs[:,2]==1)])
		# logs[:,5]	# time_taken data, filtered with the same slicing as before
		plt.subplot(2,4,3)
		plt.plot(steps_x, p_solve_y, label=_label)
		plt.subplot(2,4,4)
		plt.plot(steps_x, p_solve_y, label=_label)
		plt.subplot(2,4,7)
		plt.plot(steps_x, p_solve_y, label=_label)
		plt.subplot(2,4,8)
		p_solve_y = 1 - p_solve_y
		plt.plot(steps_x, p_solve_y, label=_label)

	for i in range(1, 9):
		plt.subplot(2,4,i)
		plt.legend()

	if _make_file is None:
		plt.show()
	else:
		plt.savefig(_make_file)
	plt.close('all')
#################################################################################
#						Class of SAT Solution Searchers							#
#################################################################################
class SATsolvers:
	def __init__(self, **kwargs):
		# Initialising all the variables required by the class
		# Attributes not available to be edited outside the class's scope
		(self.__cnf,	self.__no_v,	self.__no_c,	self.__cnf_arr,
		self.__ωp,		self.__tabu,	self.__nIters,
		self.__bIters,	self.__nFlips,	self.__type)		= [None] * 10
		self.__i, self.__f, self.__max_bI, self.__comp_ctr	= (0, 0, 0, 0)

		self.__call__(**kwargs)

	########################## User callable functions ##########################
	def __call__(self, **kwargs):
		self.__nIters	= kwargs.get('_nIters')
		self.__nFlips	= kwargs.get('_nFlips')
		self.__bIters	= kwargs.get('_bIters')
		self.__ωp		= kwargs.get('_wp')
		_tl 			= kwargs.get('_tl')
		self.__tabu		= deque(maxlen=_tl) if isinstance(_tl, int) else None
		# Creates a FIFO list of input _tl maximum size

		if '_fName' in kwargs.keys():
			(self.__cnf,	self.__no_v,	self.__no_c,
			self.__cnf_arr)	= self.__cnf_file_parser(kwargs['_fName'])

		minimum_requirement = (self.__cnf and self.__nIters
								and self.__bIters and self.__nFlips)
		# 					None	: if either of the variables is None
		searchers = {
			'GSAT'	   	: (minimum_requirement),
			'GWSAT'	   	: (minimum_requirement and self.__ωp),
			'GSAT+Tabu'	: (minimum_requirement and self.__tabu),
			'GWalkSAT' 	: (minimum_requirement and self.__ωp),
		}
		# Minimum requirements for respective search algorithms

		if '_search' in kwargs.keys():
			v = kwargs.get('_search')	# Gets the algorithm
			if v in searchers.keys():	# Checks if the algorithm is available
				self.__type = v
				if searchers[v] is not None: # if minimum requirements are met
					result = self.__local_search()
					# Performs search and stores result
					(self.__ωp,		self.__tabu,	self.__nIters,
					self.__bIters,	self.__nFlips,	self.__type) = [None] * 6
					# Resets the related variables to __init__ state
					return result
				else:
					raise Exception('Missing input variables')
			else:
				raise Exception(f'{v} Unknown Searcher')

	def grid_compare(self, _rounds, _fName, param_grid, default_kwargs={}):
		'''
		Generates csv reports of empirical results by changing one parameter
		at a time. Please check Appendix for sample code to call this function
		'''
		(self.__cnf,	self.__no_v,	self.__no_c,
		self.__cnf_arr)	= self.__cnf_file_parser(_fName)
		# Reads and parses the file for the Grid Comparison

		if not isinstance(default_kwargs, dict):
			raise Exception('default_kwargs should be a dict, review Appendix')

		if '_fName' in default_kwargs:
			del default_kwargs['_fName']

		if not isinstance(param_grid, dict):
			raise Exception('param_grid should be a dict, review Appendix')

		if '_search' in param_grid:
			_searcher = param_grid['_search']
		elif '_search' in default_kwargs:
			_searcher = default_kwargs['_search']
			del default_kwargs['_search']
		else:
			raise Exception('param_grid should have "_search" clause')
		# Gets the algorithm from either param_grid or default_kwargs

		out_file_name = f"Results_grid_compare_" + _searcher
		# eg.: Results_grid_compare_GSAT

		for item in tqdm(param_grid, desc='Parameters'):
			chk1 = re.search('_range$', item)
			# Checks if the string ends with _range
			chk2 = re.search(r'^_nFlips|^_bIters|^_nIters|^_tl|^_wp', item)
			# Checks if the string starts with names of the parameters
			if ((not chk1) and (not chk2)):
				continue
			kwargs = default_kwargs.copy()
			# Resets the kwargs (parameters for current run) to the
			# default_kwargs (default parameters)
			kwargs['_search'] = _searcher
			# sets the search algorithm
			k = item.replace('_range', '')
			empirical_log = []
			# empty list to collect empirical logs
			for v in tqdm(param_grid[item], leave=False, desc=k.ljust(10, ' ')):
				logs = []
				# empty list to collect all logs per parameter
				random.seed(myseed*v)
				# Sets random set for trackability
				for r in trange(_rounds, leave=False, desc='Iterations'):
					kwargs[k] = v

					start 			= time.time()
					result			= self.__call__(**kwargs) # Runs search
					end				= time.time()

					##################### Data Collection #######################
					time_taken		= end - start
					restarts_taken	= self.__i + 1
					flips_taken		= self.__f + 1
					bad_iter		= self.__max_bI
					sol_reached		= bool(result)
					total_runs		= restarts_taken * flips_taken
					comparison_runs	= self.__comp_ctr

					logs.append([r+1,
								sol_reached,
								restarts_taken,
								flips_taken,
								bad_iter,
								total_runs,
								comparison_runs,
								time_taken])
					#############################################################
					self.__comp_ctr = 0 # Resets counter
					self.__max_bI 	= 0 # Resets counter
				#################### Aggregation of results  ####################
				logs = np.array(logs)
				logs = logs[logs[:,1]==1]
				# Filters to only solved data
				logs_agg = [
					f'{k} = {v}',
					f'success = {logs.shape[0]},max,min,mean,median,std',
					'Restarts,',
					'Flips,',
					'BadIters,',
					'TotalIters,',
					'CalcCalls,',
					'Runtime,',
				] # List structured so as to fill a csv file
				# Sample output
				# _nFlips = 10
				# success = 65	max		min		mean		median		std
				# Restarts
				# Flips
				# BadIters
				# TotalIters
				# CalcCalls
				# Runtime


				if logs.shape[0] != 0:
					for i in range(2, 8):
						elem = [
							np.max(logs[:,i]),
	 						np.min(logs[:,i]),
							np.mean(logs[:,i]),
							np.median(logs[:,i]),
							np.std(logs[:,i])
						]
						conc_string = ','.join(map(str,elem))
						logs_agg[i] += conc_string
						# eg: 'Restarts,10,1,5.107692308,5,2.904230744'

				empirical_log.append('\n'.join(logs_agg))
				# Appends the whole 1D list into empirical_log
				#################################################################

			filename = out_file_name + f"_{k}.csv"
			data = '\n'.join(empirical_log)
			# Joins the elements of the empirical_log with a new line
			with open(filename, 'w') as f:
				f.write(data+'\n')

	def compare_searchers(self, _fName, _rounds, _searchers=[], **kwargs):
		'''
		Generates csv files of the comparison data generated here.
		Plots are drawn with a helper function plot_rtd()
		Please check Appendix for sample code to call this function.
		'''
		(self.__cnf,	self.__no_v,	self.__no_c,
		self.__cnf_arr)	= self.__cnf_file_parser(_fName)
		# Reads and parses the file for the RTD comparisons
		if '_search' in kwargs:
			del kwargs['_search']
		# Removes search algorithm from kwargs because it is to be set here

		kwargs['_nIters'] = 1 # 1 for no restarts allowed
		kwargs['_bIters'] = 0 # 0 for unlimited bad iterations
		# Setting minimum requirement variables

		if _searchers == []:
			_searchers = ['GSAT', 'GWSAT', 'GSAT+Tabu', 'GWalkSAT']
			# If the _searchers is not Initiated, then all the
			# four algorithms will be compared.

		logs = []
		for r in trange(_rounds):
			for searcher in _searchers:
				random.seed(myseed*r)
				# Setting seed per round for trackability
				kwargs['_search'] = searcher

				start 			= time.time() # Get starting time
				result			= self.__call__(**kwargs)
				end				= time.time() # Get ending time

				time_taken		= end - start # Calculate time difference
				flips_taken		= self.__f + 1 # Total flips taken to solve
				sol_reached		= bool(result) # False: if result is []
				comparison_runs	= self.__comp_ctr
				# Calls made to __comparator() function

				logs.append([r+1,
							_searchers.index(searcher),
							sol_reached,
							flips_taken,
							comparison_runs,
							time_taken])
				self.__comp_ctr = 0
				# Reset counter

		logs = np.array(logs)
		out_fName = 'Results_compare_searchers'
		np.savetxt(f'{out_fName}.csv', logs, delimiter=',', fmt='%s',
					header='   '.join([f'{i}: {_searchers[i]}'
											for i in range(len(_searchers))]) +
					'\nround,searcher,sol_reached,' +
					'flips_taken,comparison_runs,time_taken')
		# Eg Output:
		# # 0: GSAT+Tabu   1: GWalkSAT
		# # round searcher sol_reached flips_taken comparison_runs time_taken
		#	1		0			1			40			640			0.172540426
		#	1		1			1			8			36			0.007979393
		#	2		0			1			10			160			0.037935495

		print(f'Comparison results saved to "{out_fName}.csv"')

		plot_rtd(logs, _use_legend=_searchers, _make_file=f'{out_fName}.png')
	############################# Internal functions ############################
	def __cnf_file_parser(self, file_name):
		with open(file_name, 'r') as f:
			txt = f.read().strip()
		# Open file --> Read file --> Strip spaces --> Store --> Close file

		while '  ' in txt:
			txt = txt.replace('  ', ' ')
		# Replace redundant double spaces (if any) to single spaces

		first_split	= [
			line.strip()
			for line in txt.split(sep='\n')
			if line.strip()[0].lower() not in ['c', '%', '']
		]
		# Split lines of text which does not start with c, %, or blank, eg: -
		#
		# c    mixed sat? no
		# c    clause length = 3
		# p cnf 20  91
		#  6 -3 13 0			--> ['p cnf 20  91', '6 -3 13', '5 9 -14', '0']
		# 5 9 -14 0
		# %
		# 0

		if first_split[0][0] == 'p':
			v, c = map(int, first_split[0].split(sep=' ')[2:])
			# 'p cnf 20 91' -> ['p', 'cnf, '20', '91'] -> [20, 91]
		else:
			raise Exception('Not a SAT instance file')

		clauses = ' '.join(first_split[1:]).strip(' 0').split(sep=' 0 ')
				# join the list with 		strip the	split the whole by
				# spaces in between 		ending 0	' 0 '
			# eg: ..8 -14 0 13 -11..					..8 -14'], ['13 -11..
		clauses = [set(map(int, c_str.split(sep=' '))) for c_str in clauses]
			# eg: ...['-19 -12 -15'], ['-5 -15 -18'], ['5 9 -14'],... -->
			#	 [...{-19, -12, -15}, {-5, -15, -18}, {5, 9, -14},...]

		return clauses, v, c, np.array(list(map(list, clauses)))

	def __sol_encoder(self, sol):
		'''
		Takes in bit coded sol returns variable encoded sol
		sol format should be a index-wise bit coded list
		Index = variable name, Value = 1 or 0 (meaning True or False)
		sol[0] should be None, since there is no variable named 0
		eg: [None, 0, 1, 0] --> [-1, 2, -3]
		'''
		if (sol[0] == None) and not bool(set(sol[1:]) - {0, 1}):
		#	checks Index 0	checks if rest are only 0 or 1
			return {(-1 * i) if sol[i] == 0 else i for i in range(1, len(sol))}
				#	-1 x Index if value is 0
				#				otherwise the Index itself
		else:
			raise Exception(f'{sol} is not an index-wise bit coded list')

	def __comparator(self, cnf_set, sol_set, encode=False, net_gain=False, sat_unsat_bool=False):
		'''
		Takes in the cnf clauses and solution and return True if the solution
		satisfies the all the cnf clauses.
		solution format: [-1, -5, ..... 19]

		net_gain = True, returns % of variables satisfied
		sat_unsat_bool = True, returns an array True(sat)/False(unsat) clause
		encode parameter is to support solution in binary format
		[None, 0, 0, 1 ... 1, 0]
		'''
		if encode:
			sol_set = self.__sol_encoder(sol_set)
			# gets the bit-coded sol_set encoded to varible-state list

		self.__comp_ctr += 1
		if net_gain or sat_unsat_bool:
			sat_unsat_arr = np.array(
						[
							True if clause.intersection(sol_set) else False
							# If the clause has any of the variables in
							# sol_set, then add True into the list
							for clause in cnf_set
						]
					)
			if sat_unsat_bool:
				return sat_unsat_arr
			else:
				return sum(sat_unsat_arr) * 100 / self.__no_c
						# Divide by length of cnf_set, x 100 for %
						# Output: 100.0 means sol_set is a solution

		for clause in cnf_set:
			if not clause.intersection(sol_set):
				# If the clause does not have any of the variable in sol_set
				return False
		return True

	def __select_variable_GSAT(self, A, _collection=None):
		net_gain = []
		if _collection is None:
			_collection = range(1, len(A))
		# Option for other select_variable functions to resuse this
		# for GSAT/GWSAT the loop would range from 1 to n number of variables
		# for GSAT_Tabu the collection is a set of varibles except the tabu
		# for GWalkSAT, variables present in the unsolved clause

		for i in _collection:
			sol = A.copy()
			sol[i] = int(not(sol[i]))
					# Flips variable at i-th position
			gain = self.__comparator(self.__cnf, sol, encode=True, net_gain=True)
			# Gets gain for the solution with the varible at i-th position flipped
			net_gain.append([i, gain])
			# Stores for comparison later

		net_gain = np.array(net_gain)
		# Converts to array for easy slicing/flitering
		max_gain = net_gain[:,1].max()
		# Gets the max gain in the array
		max_gain_indices = net_gain[net_gain[:,1] == max_gain][:,0]
		# net_gain[:,1] == max_gain	# Filters indices with max gain
		# Case: Multiple indices with same gain
		# Column 0 is same as index for GSAT/GWSAT
		# but different for Tabu and GWalkSAT
		choice = int(random.choice(max_gain_indices))
		# Randomly chooses one
		return choice

	def __select_variable_GWSAT(self, A):
		if random.random() < self.__ωp:
			# Random walk step
			sat_unsat_bool = self.__comparator(self.__cnf, A,
											   encode=True, sat_unsat_bool=True)
			# Array of [True, False, ....True] for [Sat, unsat, ....Sat] clauses
			unsat_clauses = self.__cnf_arr[sat_unsat_bool==False]
			# Filters only the unsat clauses
			unsat_vars = np.unique(unsat_clauses.flatten())
			# Flattens to an array of unique variables involved in unsat clauses
			choice = random.choice(unsat_vars)
			# Randomly chooses variables
			return choice
		else:
			return self.__select_variable_GSAT(A)

	def __select_variable_GSAT_Tabu(self, A):
		sol_indices_to_check = set(range(1,len(A))) - set(self.__tabu)
		# All the indices (variables) except the ones tabu
		choice = self.__select_variable_GSAT(A, sol_indices_to_check)
		self.__tabu.append(choice)
		# Insers the chosen variable on the right and pops a variable from the
		# left to maintain the length of TABU to _tl
		return choice

	def __select_variable_GWalkSAT(self, A):
		sat_unsat_bool = self.__comparator(self.__cnf, A,
											   encode=True, sat_unsat_bool=True)
		# Array of [True, False, ....True] for [Sat, unsat, ....Sat] clauses
		unsat_clauses = self.__cnf_arr[sat_unsat_bool==False]
		# Filters only the unsat clauses
		unsat_clause = random.choice(unsat_clauses)
		# Randomly chooses an unsat clause
		if random.random() < self.__ωp:
			# Random walk step
			unsat_var = abs(random.choice(unsat_clause))
			# Randomly choose a variable from unsat clause
			# abs() because NOT operator is not part of the varible name
			# eg: In -17, 17 is the variable name
			return unsat_var
		else:
			sol_indices_to_check = list(map(abs, unsat_clause))
			# To get gains for all the variables in the clause
			return self.__select_variable_GSAT(A, sol_indices_to_check)

	def __select_variable_dispatcher(self, A, type):
		'''
		Calling the respective select_variable function, based on the algorithm
		'''
		return {
			'GSAT'		: self.__select_variable_GSAT,
			'GWSAT'		: self.__select_variable_GWSAT,
			'GSAT+Tabu'	: self.__select_variable_GSAT_Tabu,
			'GWalkSAT'	: self.__select_variable_GWalkSAT,
		}.get(type)(A)
		# Calls the function with parameter A, returns None if the algorithm
		# doesn't exist

	def __local_search(self):
		'''General local search function to be used for all the algorithms'''
		for self.__i in range(self.__nIters):
			A = [None] + [random.choice([0,1]) for _v in range(self.__no_v)]
			# Initialising a random solution. Indexes define the variable
			# Values define if it is True/False, zeroth position is None because
			# there is no variable called 0
			# eg.: [None, 0, 1, 0, 0, ....]
			last_gain = 0		# Variable dedicated to take care of the
			bad_iter_ctr = 0	# stagnation parameter _bIters
			for self.__f in range(self.__nFlips):
				gain = self.__comparator(self.__cnf, A,
										encode=True, net_gain=True)
				# Gets and stores the gain proposed solution
				if gain == 100.0:
					# the gain is at 100%, returns the found solution
					sol = self.__sol_encoder(A)
					# Converts binary True/False to variables as numbers
					# eg: [None, 0, 1, 0, 0, ...] -> [-1, 2, -3, -4, ...]
					sol = sorted(sol, key=abs)
					# eg: [-19, -17, ... 5, 7, 9] -> [-1, ... 7, 9....-17, -19]
					return sol
				if (self.__bIters is not None) and (self.__bIters > 0):
					bad_iter_ctr  = bad_iter_ctr + 1 if gain == last_gain else 0
					# If the last found gain is same as the current found gain
					# bad_iter_ctr is incremented, otherwise set as 0
					self.__max_bI = (bad_iter_ctr
									if bad_iter_ctr > self.__max_bI
									else self.__max_bI)
					# Max bad iter while the whole solution finding, used for
					# Empirical results collation
					if bad_iter_ctr == self.__bIters:
						# If bad counter reaches allowed limit, break the loop
						break
					last_gain = gain
				x = self.__select_variable_dispatcher(A, self.__type)
				# Gets the best variable to flip
				A[x] = int(not A[x])
				# Flips the varible
		return []
		# Returns empty list if solution is not found
#################################################################################
#								sanity checker									#
#################################################################################
def sanity(inst_file, alg, nRuns, nBadIters, nRestarts, wp, tl):
	try:
		if not os.path.isfile(inst_file):
			# checks if the instance path is a file
			raise Exception()

		alg = {
			'GSAT'		: 'GSAT',
			'GWSAT' 	: 'GWSAT',
			'GSATTABU'	: 'GSAT+Tabu',
			'GTABU'		: 'GSAT+Tabu',
			'GWALKSAT'	: 'GWalkSAT'
		}.get(alg.upper())
		# Changes the inputs to the expected inputs by the class

		if alg not in ['GSAT','GWSAT','GSAT+Tabu','GWalkSAT']:
			# Checks if the algorithm exists in class's capabilities
			raise Exception()

		nRuns		= int(nRuns)		# Casting as integer
		nBadIters	= int(nBadIters)	# Casting as integer
		nRestarts	= int(nRestarts)	# Casting as integer
		wp			= float(wp)			# Casting as float
		if not 0 <= wp <= 1:			# wp must be a value between 0 and 1
			raise Exception()
		tl			= int(tl)			# Casting as integer

		return inst_file, alg, nRuns, nBadIters, nRestarts, wp, tl
		# Returns the converted inputs
	except:
		print('Problem in command line arguments')
		exit(0)
#################################################################################
#									__main__									#
#################################################################################
if __name__ == '__main__':
	if len(sys.argv) < 8:
		print('Frame the command line arguments correctly again')
		exit(0)
	inst_file, alg, nRuns, nBadIters, nRestarts, wp, tl = sanity(*sys.argv[1:])
	myseed = 195638
	start = time.time()
	sat = SATsolvers()
	sol = sat(
			_fName=inst_file,
			_nIters=nRestarts,
			_bIters=nBadIters,
			_nFlips=nRuns,
			_wp=wp,
			_tl=tl,
			_search=alg
			)
	end = time.time()
	print(sol,
		'Execution time: '+str(round((end - start), 4))+' seconds', sep='\n\n')
#################################################################################
#								Appendix: Code Samples							#
#################################################################################
'''
Explains the code samples for user callable functions in the class.
'''
if False:
	fName = [
		r'./uf20_91/inst/uf20-0385.cnf',
		r'./uf20_91/inst/uf20-0876.cnf',
		r'./uf20_91/inst/uf20-0961.cnf'
	]
	cwd = os.getcwd()
	fName = [cwd + f[1:] for f in fName]
	myseed = 195638
	sat = SATsolvers()
	# To find solution to an instance
	sat(_fName=fName[0], 	# instance path
		_nIters=10,			# allowed restarts
		_bIters=10,			# allowed number of non-improving iterations
		_nFlips=100,		# allowed number of iterations per restart
		_wp=0.20,			# walk probability
		_tl=5,				# tabu length
		_search='GSAT'		# to trigger the search algorithm
		)

	# To polt the RTD graphs, if the results are already evaluated
	plot_rtd(logs=None,
			_use_file='Results_compare_searchers.csv',
			_use_legend=True)

	# To perform RTD plots
	sat.compare_searchers(
		_rounds=10000,
		_fName=fName[2],
		_searchers=['GWSAT', 'GSAT+Tabu'],
		_nFlips=1000,
		_wp=0.20,
		_tl=5
	)

	# Experimental setups to produce Empirical Results into csv files
	sat.grid_compare(
		_rounds=100,
		_fName=fName[0],
		param_grid={
				'_search'		: 'GSAT',
				'_nIters_range' : range(10,100,10),
				'_nFlips_range'	: range(10,100,10),
				'_bIters_range'	: range(5,30,5),
				'_wp_range'		: np.linspace(0.1,0.5,5),
				'_tl_range'		: range(1,6,1),
			},
		default_kwargs={
				'_search'	: 'GSAT',
				'_nIters'	: 10,
				'_bIters'	: 20,
				'_nFlips'	: 50,
				'_wp'		: 0.1,
				'_tl'		: 5
			}
		)
