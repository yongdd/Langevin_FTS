import os
import time
import numpy as np
from scipy.io import savemat
import scft

# OpenMP environment variables
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"  # 0, 1
os.environ["OMP_NUM_THREADS"] = "2"  # 1 ~ 4

# Major Simulation params
f = 0.5         # A-fraction of major BCP chain, f
eps = 1.0       # a_A/a_B, conformational asymmetry

params = {
    # "platform":"cpu-mkl",           # choose platform among [cuda, cpu-mkl]
    
    "nx":[32],          # Simulation grid numbers
    "lx":[4.36],        # Simulation box size as a_Ref * N_Ref^(1/2) unit,
                        # where "a_Ref" is reference statistical segment length
                        # and "N_Ref" is the number of segments of reference linear homopolymer chain.

    "box_is_altering":True,     # Find box size that minimizes the free energy during saddle point iteration.
    "chain_model":"discrete",   # "discrete" or "continuous" chain model
    "ds":1/90,                  # Contour step interval, which is equal to 1/N_Ref.

    "segment_lengths":{         # Relative statistical segment length compared to "a_Ref.
        "A":np.sqrt(eps*eps/(eps*eps*f + (1-f))), 
        "B":np.sqrt(    1.0/(eps*eps*f + (1-f))), },

    "chi_n": [["A", "B", 13.27]],   # Interaction parameter, Flory-Huggins params * N_Ref

    "distinct_polymers":[{      # Distinct Polymers
        "volume_fraction":1.0,  # volume fraction of polymer chain
        "blocks":[              # AB diBlock Copolymer
            {"type":"A", "length":f, }, # A-block
            {"type":"B", "length":1-f}, # B-block
        ],},],

    "am":{
        "max_hist":20,          # Maximum number of history
        "start_error":1e-2,     # When switch to AM from simple mixing
        "mix_min":0.1,          # Minimum mixing rate of simple mixing
        "mix_init":0.1,         # Initial mixing rate of simple mixing
    },

    "max_iter":2000,     # The maximum relaxation iterations
    "tolerance":1e-8     # Terminate iteration if the self-consistency error is less than tolerance
}

# Set initial fields
w_A = np.zeros(list(params["nx"]), dtype=np.float64)
w_B = np.zeros(list(params["nx"]), dtype=np.float64)
print("w_A and w_B are initialized to lamellar phase.")
for i in range(0,params["nx"][0]):
    w_A[i] =  np.cos(3*2*np.pi*i/params["nx"][0])
    w_B[i] = -np.cos(3*2*np.pi*i/params["nx"][0])

# Initialize calculation
calculation = scft.SCFT(params=params)

# Set a timer
time_start = time.time()

# Run
calculation.run(initial_fields={"A": w_A, "B": w_B})

# Estimate execution time
time_duration = time.time() - time_start
print("total time: %f " % time_duration)

# Save final results
phi = calculation.get_concentrations()
w = calculation.get_fields()

mdic = {"params":params, "dim":len(params["nx"]), "nx":params["nx"], "lx":params["lx"], "ds":params["ds"],
        "f":f, "chi_n":params["chi_n"], "epsilon":eps, "chain_model":params["chain_model"],
        "w_a":w_A, "w_b":w_B, "phi_a":phi["A"], "phi_b":phi["B"]}
savemat("fields.mat", mdic)

# Recording first a few iteration results for debugging and refactoring
    #    1    2.220E-16  [ 1.0481009E+00  ]    -0.009300846   2.2442714E-01  [  4.3600000 ]
    #    2    0.000E+00  [ 1.0505473E+00  ]    -0.009752572   2.2487508E-01  [  4.3597350 ]
    #    3    0.000E+00  [ 1.0530913E+00  ]    -0.010221266   2.2532266E-01  [  4.3594606 ]
    #    4   -6.661E-16  [ 1.0557620E+00  ]    -0.010712001   2.2576538E-01  [  4.3591738 ]
    #    5    0.000E+00  [ 1.0585655E+00  ]    -0.011225570   2.2619631E-01  [  4.3588742 ]
