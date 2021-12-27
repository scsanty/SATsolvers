# SATsolvers
Solving SAT problems using GSAT, GWSAT, GSAT using Tabu, and GWalkSAT

Libraries: -
1) numpy is required

2) matplotlib is used for plots (optional)

3) tqdm is (optional)

python <Chowdhury_Sayantan_R00195638_MH2.py> <inst file> <alg> <nRuns> <nBadIters> <nRestarts> <wp> <tl>

where

- inst file: the full file path for the instance

- alg: one of f gsat, gwsat, gtabu, gwalksat

- nRuns: the number of runs of the algorithm to perform

- nBaditers: is the number of successive iterations without an improving move before restarting

- nRestarts: is the number of restarts

- wp: is the random walk probability

- tl: the length of the tabu list

REM SAMPLES

ECHO GSAT

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0385.cnf gsat 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0876.cnf gsat 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0961.cnf gsat 10 20 50 0.1 5

ECHO GWSAT

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0385.cnf gwsat 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0876.cnf gwsat 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0961.cnf gwsat 10 20 50 0.1 5

ECHO GSAT+Tabu

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0385.cnf gtabu 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0876.cnf gtabu 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0961.cnf gtabu 10 20 50 0.1 5

ECHO GWalkSAT

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0385.cnf gwalksat 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0876.cnf gwalksat 10 20 50 0.1 5

python Chowdhury_Sayantan_R00195638_MH2.py uf20_91/inst/uf20-0961.cnf gwalksat 10 20 50 0.1 5

N.B. 

1) Outputs of the experiments are stored in Output folder

2) uf20_91 folders contains the instances

3) helper-plots.py was used to plot the graphs shown in the Report, and is not connected to the main code in any manner

############################# Appendix: Code Samples ############################

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

sat(_fName=fName[0], # instance path

_nIters=10,  # allowed restarts

_bIters=10,  # allowed number of non-improving iterations

_nFlips=100,  # allowed number of iterations per restart

_wp=0.20,  # walk probability

_tl=5,  # tabu length

_search='GSAT'  # to trigger the search algorithm

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

'_search'  : 'GSAT',

'_nIters_range' : range(10,100,10),

'_nFlips_range'  : range(10,100,10),

'_bIters_range'  : range(5,30,5),

'_wp_range'  : np.linspace(0.1,0.5,5),

'_tl_range'  : range(1,6,1),

},

default_kwargs={

'_search'  : 'GSAT',

'_nIters'  : 10,

'_bIters'  : 20,

'_nFlips'  : 50,

'_wp'  : 0.1,

'_tl'  : 5

}

)
