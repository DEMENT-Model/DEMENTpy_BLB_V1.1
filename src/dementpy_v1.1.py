"""
-------------------------------------------------------------------------------
      DEMENTpy--Decomposition Model of Enzymatic Traits in Python,v1.1
                              Bin Wang, Ph.D.
              Department of Ecology and Evolutionary Biology
                       University of California Irvine
                  Emails: wbwenwu@gmail.com or bwang7@uci.edu
                          Twitter: @bio_atmosphere
                          
                  Initialization processes added 2022-03-11
                          Brittni Bertolet, Ph.D.
-------------------------------------------------------------------------------
"""
import os
import sys
import pandas as pd
import numpy  as np
import pickle

from grid import Grid
from output import Output
from utility import export

from substrate import Substrate
from monomer   import Monomer
from enzyme    import Enzyme
from microbe   import Microbe
from utility   import expand

#...Define function for running DEMENT with either a new grid or a previously initialized grid 
def main():
    
    
    print("""
    ---------------------------------------------------------------------------
         DEMENTpy (DEcomposition Model of Enzymatic Traits in Python)
                               Version 1.0
               Department of Ecology and Evolutionary Biology
                     University of California Irvine
    ---------------------------------------------------------------------------       
    """)
    
     
    #...Obtain the command line arguments
    input_folder  = sys.argv[1]   # input folder name
    output_folder = sys.argv[2]   # output folder name
    outname       = sys.argv[3]   # output file name
    
    #...Obtain the command line arguments for initializing the grid
    grid_init     = bool(sys.argv[4])    # True/False for grid initialization
    grid          = sys.argv[5]    # grid file name to read in
    
    #...Set up the working directory
    os.chdir('../'+input_folder)

    #...grow a seed of random number generator
    np.random.seed(int(outname[:-4]))

    #...a few system constants
    runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
    pulse      = int(runtime.loc['pulse',1])         # number of pulses
    cycle      = int(runtime.loc['end_time',1])      # number of time steps in each pulse
    interval   = int(runtime.loc['interval',1])      # interval of time step to record outputs
    mic_reinit = True    # indicate reinitialization of microbial community
    
    #...Set up initialization
    # Load all input files
    parameters      = pd.read_csv('parameters.csv',         header=None, index_col=0).astype('float32')   # parameters
    substrates_init = pd.read_csv('initial_substrates.csv', header=0,    index_col=0).astype('float32')   # initial substrates
    sub_mon_input   = pd.read_csv('sub_mon_inputs.csv',     header=0,    index_col=0).astype('float32')   # inputs of substrates and monomers
    Ea_input        = pd.read_csv("enzyme_ea.csv",          header=0,    index_col=0).astype('float32')   # enzyme activation energy
    climate         = pd.read_csv('climate.csv',            header=0,    index_col=0)                     # climate forcings

    # daily temperature and water potential
    daily_temp = climate['Temp'].astype('float32')  # temperaure series
    daily_psi  = climate['Psi'].astype('float32')   # water potential series

    if grid_init == True: # Initialize a new grid 
        #...an instance of Substrate class 
        Substrates = Substrate(runtime,parameters,substrates_init)
        #...substrate initial pool size
        substrates_initial_pool = Substrates.Substrates_start
        #...substrate input rate
        substrates_input_rate = Substrates.substrate_input(sub_mon_input)
        #...substrates-produced monomers
        substrates_produced_monomers = Substrates.substrate_produced_monomer()
        #...substrates degradation required enzymes
        substrates_req_enzyme = Substrates.substrate_degradation_enzyme()

        #...an instance of Monomer class
        Monomers = Monomer(runtime,parameters)
        #...monomers initial pool size
        monomers_initial_pool = Monomers.monomer_initialization(substrates_initial_pool)
        #...initial monomer ratios
        monomer_ratio_inital = Monomers.monomer_ratios(monomers_initial_pool)
        #...monomers input rate
        monomers_input_rate = Monomers.monomer_input_rate(sub_mon_input)
        #...monomers uptake required enzymes
        monomers_uptake_reqenzyme = Monomers.monomer_uptake_reqenzyme()

        #...an instance of Enzyme class
        Enzymes = Enzyme(runtime,parameters,substrates_initial_pool.index)
        #...enzyme initial pool size:0
        enzymes_initial_pool = Enzymes.enzyme_pool_initialization()
        #...enzyme attributes
        enzymes_attributes = Enzymes.enzyme_attributes()
        #...enzymes of substrate degradation Ea    
        enzymes_Ea = Enzymes.enzyme_Ea(Ea_input)
        #...monomers uptake enzyme Ea
        enzymes_uptake_Ea = Enzymes.enzyme_uptake_Ea()
        #...enzymes of substrate degradation Vmax
        enzymes_Vmax,enzymes_Vmax_T = Enzymes.enzyme_Vmax(substrates_req_enzyme)
        #...monomers uptake enzyme Vmax
        enzymes_uptake_Vmax= Enzymes.enzyme_uptake_Vmax(monomers_uptake_reqenzyme)
        #...enzymes of substrate degradation Km
        enzymes_Km = Enzymes.enzyme_Km(enzymes_Vmax)
        #...monomers uptake enzyme Km
        enzymes_uptake_Km = Enzymes.enzyme_uptake_Km(enzymes_uptake_Vmax)

    else: # load in a previously initialized grid 
        with open(str(grid), "rb") as fh:
            grid_initialization = pickle.load(fh) # write a function that initializes one grid

        #...substrate initial pool size
        substrates_initial_pool = grid_initialization['Substrates']
        #...substrate input rate
        substrates_input_rate = grid_initialization['SubInput']
        #...substrates-produced monomers
        substrates_produced_monomers = grid_initialization['MonomersProduced']
        #...substrates degradation required enzymes
        substrates_req_enzyme = grid_initialization['ReqEnz']

        #...monomers initial pool size
        monomers_initial_pool = grid_initialization['Monomers']
        #...initial monomer ratios
        monomer_ratio_inital = grid_initialization['Monomer_ratio']
        #...monomers input rate
        monomers_input_rate = grid_initialization['MonInput']
        #...monomers uptake required enzymes
        monomers_uptake_reqenzyme = grid_initialization['Uptake_ReqEnz']

        #...enzyme initial pool size:0
        enzymes_initial_pool = grid_initialization['Enzymes']
        #...enzyme attributes
        enzymes_attributes = grid_initialization['EnzAttrib']
        #...enzymes of substrate degradation Ea    
        enzymes_Ea = grid_initialization['Ea']
        #...monomers uptake enzyme Ea
        enzymes_uptake_Ea = grid_initialization['Uptake_Ea']
        #...enzymes of substrate degradation Vmax
        enzymes_Vmax_T = grid_initialization['Vmax0']
        #...monomers uptake enzyme Vmax
        enzymes_uptake_Vmax= grid_initialization['Uptake_Vmax0']
        #...enzymes of substrate degradation Km
        enzymes_Km = grid_initialization['Km0']
        #...monomers uptake enzyme Km
        enzymes_uptake_Km = grid_initialization['Uptake_Km0']

    # Initialize the microbial community 
    #...an instance of Microbe class
    Microbes = Microbe(runtime,parameters)
    #...Microbial community initialization#...note microbial_community is a tuple
    microbial_community = Microbes.microbial_community_initialization()
    #...Microbial minimum ratios
    microbial_min_ratios = Microbes.minimum_cell_quota()
    #...Microbial enzyme genes
    microbial_enzyme_gene = Microbes.microbe_enzyme_gene()
    #...Microbial osmolyte genes
    microbial_osmolyte_gene = Microbes.microbe_osmolyte_gene()
    #...Microbial uptake genes
    microbial_uptake_gene = Microbes.microbe_uptake_gene(substrates_req_enzyme,microbial_enzyme_gene,substrates_produced_monomers)
    #...Microbial uptake cost
    microbial_uptake_cost = Microbes.microbe_uptake_cost(microbial_uptake_gene)
    #...Microbial enzyme production rate
    microbial_enzyme_prod_rate = Microbes.microbe_enzproduction_rate(microbial_enzyme_gene,enzymes_attributes)
    #...Microbial osmolyte productoin rate
    microbial_osmolyte_prod_rate = Microbes.microbe_osmoproduction_rate(microbial_osmolyte_gene)
    #...Microbial drought tolerance
    microbial_drought_tol = Microbes.microbe_drought_tol(microbial_osmolyte_prod_rate[2],microbial_osmolyte_prod_rate[3])
    #...Microbial mortality
    microbial_mortality = Microbes.microbe_mortality(microbial_community[2])

    # Put it in a dictionary
    microbe_initializtion = {"microbial_community":     microbial_community,
                             "microbial_min_ratios":    microbial_min_ratios,
                             "microbial_enzyme_gene":   microbial_enzyme_gene,
                             "microbial_osmolyte_gene": microbial_osmolyte_gene,
                             "microbial_uptake_gene":   microbial_uptake_gene,
                             "microbial_uptake_cost":   microbial_uptake_cost,
                             "microbial_enzyme_prod_rate":   microbial_enzyme_prod_rate,
                             "microbial_osmolyte_prod_rate": microbial_osmolyte_prod_rate,
                             "microbial_drought_tol":   microbial_drought_tol,
                             "microbial_mortality":     microbial_mortality
    }

    #...Dump all initialized data into a dictionary; NOTE: variables with expand() put on the spatial grid
    gridsize = int(runtime.loc['gridsize',1]) # define grid size
    data_initialization = {"Substrates": expand(substrates_initial_pool,gridsize),
                    "SubInput":   expand(substrates_input_rate,gridsize),
                    "ReqEnz":           substrates_req_enzyme,
                    "MonomersProduced": substrates_produced_monomers,
                    "Monomers":     expand(monomers_initial_pool,gridsize),
                    "Monomer_ratio":expand(monomer_ratio_inital,gridsize),
                    "MonInput":     expand(monomers_input_rate,gridsize),
                    "Uptake_ReqEnz":expand(monomers_uptake_reqenzyme,gridsize),
                    "Enzymes":      expand(enzymes_initial_pool,gridsize),
                    "Km0":          expand(enzymes_Km,gridsize),            # enzyme half-saturation constant
                    "Uptake_Km0":   expand(enzymes_uptake_Km,gridsize),     # transporter half-saturation constant
                    "Uptake_Ea":    expand(enzymes_uptake_Ea,gridsize),     # transporter acitivation energy
                    "Uptake_Vmax0": expand(enzymes_uptake_Vmax,gridsize),   # transporter reaction rate
                    "Ea":           expand(enzymes_Ea,gridsize),            # enzyme activation energy
                    "Vmax0":        expand(enzymes_Vmax_T,gridsize),        # enzyme reaction rate
                    "EnzAttrib":    enzymes_attributes,                     # enzyme stoichiometry and energy cost
                    "Microbes_pp": microbial_community[0],                  # tuple[0]: microbes preceding placement
                    "Microbes":    microbial_community[1],                  # tuple[1]: initialized spatial microbes
                    "fb":          microbial_community[2],                  # tuple[2]: fungi index
                    "Bac_density": microbial_community[3],                  # tuple[3]: bacterial density
                    "Fun_density": microbial_community[4],                  # tuple[4]: fungi density
                    "MinRatios":   expand(microbial_min_ratios,gridsize),   # microbial cell min. ratios
                    "UptakeGenes": expand(microbial_uptake_gene,gridsize),  # transporter gene distribution across taxa
                    "OsmoGenes":   expand(microbial_osmolyte_gene,gridsize),# osmolyte gene distribution across taxa
                    "EnzGenes":    expand(microbial_enzyme_gene,gridsize),  # enzyme gene distribution across taxa
                    "UptakeGenes_trait":   expand(microbial_uptake_cost[0],gridsize),        # single gene cost of transporter
                    "OsmoProdConsti_trait":expand(microbial_osmolyte_prod_rate[0],gridsize), # single gene cost of constitutive osmolyte
                    "OsmoProdInduci_trait":expand(microbial_osmolyte_prod_rate[1],gridsize), # single gene cost of inducible osmolyte
                    "EnzProdConsti_trait": expand(microbial_enzyme_prod_rate[0],gridsize),   # single gene cost of constitutive enzyme
                    "EnzProdInduci_trait": expand(microbial_enzyme_prod_rate[1],gridsize),   # single gene cost of inducible enzyme
                    "UptakeGenesCost":     expand(microbial_uptake_cost[1],gridsize),        # distribution of transporter gene cost across taxa
                    "OsmoProdConsti":      expand(microbial_osmolyte_prod_rate[2],gridsize), # distribution of consti. osmolyte gene cost across taxa
                    "OsmoProdInduci":      expand(microbial_osmolyte_prod_rate[3],gridsize), # distribution of induci. osmolyte gene cost across taxa
                    "EnzProdConstit":      expand(microbial_enzyme_prod_rate[2],gridsize),   # distribution of consti. enzyme gene cost across taxa
                    "EnzProdInduce":       expand(microbial_enzyme_prod_rate[3],gridsize),   # distribution of induci. enzyme gene cost across taxa
                    "TaxDroughtTol":       expand(microbial_drought_tol,gridsize),           # distribution of taxon-specific drought tol.
                    'basal_death_prob':  microbial_mortality[0],                # basal death probability
                    'death_rate':        microbial_mortality[1],                # change rate of death prob. agaist mositure
                    "AE_ref":            parameters.loc["CUE_ref",1],           # Reference assimilation efficiency: 0.5
                    "AE_temp":           parameters.loc["CUE_temp",1],          # AE temperature sensitivity; default: -0.016
                    'Uptake_Maint_cost': parameters.loc['Uptake_Maint_cost',1], # constant of transporter maintenence cost
                    'C_min':             parameters.loc['C_min',1],             # C threshold of cell lysis
                    'N_min':             parameters.loc['N_min',1],             # N threshold of cell lysis
                    'P_min':             parameters.loc['P_min',1],             # P threshold of cell lysis
                    'max_size_b':        parameters.loc['max_size_b',1],        # C quota threshold for bacterial cell division
                    'max_size_f':        parameters.loc['max_size_f',1],        # C quota threshold for fungal cell division
                    'wp_fc':             parameters.loc['wp_fc',1],             # threshold below which microbes start to respond to drought
                    'wp_th':             parameters.loc['wp_th',1],             # threshold below which microbes in full swing to respond to drought
                    'alpha':             parameters.loc['alpha',1],             # factor delineating curve concavity of microbial response to drought
                    'Temp': daily_temp,                                         # temperature
                    'Psi':  daily_psi                                           # water potential
                  }


    #...Prepare for output by creating an instance of the Output class
    Output_init = Output(runtime,data_initialization, microbe_initializtion)

    #...Create an instance of the Grid class
    Ecosystem = Grid(runtime,data_initialization)

    #...Run the model
    for p in range(pulse):
        
        for i in range(p*cycle, (p+1)*cycle):
        
            # substrates degradation
            Ecosystem.degradation(p,i)
        
            # monomers uptake
            Ecosystem.uptake(p,i)
        
            # microbial metabolism
            Ecosystem.metabolism(i)
        
            # microbial death
            Ecosystem.mortality(i)
        
            # microbial reproduction and dispersal
            Ecosystem.reproduction(i)
        
            # output data using the "output" method in the Output class
            if i == 0:
                Output_init.output(Ecosystem,i)  # day 1
            elif i%interval==interval-1:
                Output_init.output(Ecosystem,i)  # interval
            
            # if only 1 pusle, skip all following lines within this loop
            #if pulse == 1:
            #    continue
            
            # output microbial mass of every iteration using the "microbes_df" method in the Output class
            Output_init.microbes_abundance(Ecosystem,i)
            
            # re-initialize microbial community in each new pulse
            if i == (p+1)*cycle-1:
                Ecosystem.repopulation(Output_init,i,mic_reinit)
    
    #...export the Output_init object to the output_folder using the export() funtion in the utility module 
    os.chdir('../'+output_folder)
    export(Output_init,outname)
   
main()
