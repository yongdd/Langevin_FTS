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

chi_n = 13.27
a_ratio = 0.25
b_ratio = (1-np.sqrt(a_ratio))**2

params = {
    # "platform":"cuda",           # choose platform among [cuda, cpu-mkl]
    
    "nx":[32,32,32],          # Simulation grid numbers
    "lx":[4.36,4.36,4.36],    # Simulation box size as a_Ref * N_Ref^(1/2) unit,
                              # where "a_Ref" is reference statistical segment length
                              # and "N_Ref" is the number of segments of reference linear homopolymer chain.

    "box_is_altering":False,      # Find box size that minimizes the free energy during saddle point iteration.
    "chain_model":"continuous",   # "discrete" or "continuous" chain model
    "ds":1/100,                   # Contour step interval, which is equal to 1/N_Ref.

    "segment_lengths":{         # Relative statistical segment length compared to "a_Ref.
        "A":1.0,
        "B":1.0,
        "C":1.0,},
 
    "chi_n": {"A,B": chi_n,     # Interaction parameter, Flory-Huggins params * N_Ref
              "A,C": chi_n*(a_ratio ),
              "B,C": chi_n*(b_ratio ),},

    "distinct_polymers":[{      # Distinct Polymers
        "volume_fraction":1.0,  # volume fraction of polymer chain
        "blocks":[              # AB diBlock Copolymer
            {"type":"A", "length":0.5}, # A-block
            {"type":"B", "length":0.25}, # B-block
            {"type":"C", "length":0.25}, # C-block
        ],},],

    "optimizer":{
        "name":"am",            # Anderson Mixing
        "max_hist":20,          # Maximum number of history
        "start_error":1e-2,     # When switch to AM from simple mixing
        "mix_min":0.1,          # Minimum mixing rate of simple mixing
        "mix_init":0.1,         # Initial mixing rate of simple mixing
    },

    "max_iter":2,     # The maximum relaxation iterations
    "tolerance":1e-8     # Terminate iteration if the self-consistency error is less than tolerance
}

# Set initial fields
# w_A = np.zeros(list(params["nx"]), dtype=np.float64)
# w_B = np.zeros(list(params["nx"]), dtype=np.float64)
# w_C = np.zeros(list(params["nx"]), dtype=np.float64)
# print("w_A and w_B are initialized to lamellar phase.")
# for i in range(0,params["nx"][2]):
#     w_A[:,:,i] =  np.cos(3*2*np.pi*i/params["nx"][2])
#     w_B[:,:,i] = -np.cos(3*2*np.pi*i/params["nx"][2])
#     w_C[:,:,i] = -np.cos(3*2*np.pi*i/params["nx"][2])

w_A = np.random.normal(0.0, 5.0, params["nx"])
w_B = np.random.normal(0.0, 5.0, params["nx"])
w_C = np.random.normal(0.0, 5.0, params["nx"])

# Initialize calculation
calculation = scft.SCFT(params=params)

# Set a timer
time_start = time.time()

# Run
calculation.run(initial_fields={"A": w_A, "B": w_B, "C": w_C})

# Estimate execution time
time_duration = time.time() - time_start
print("total time: %f " % time_duration)

# Save final results
calculation.save_results("fields.mat")

# Recording first a few iteration results for debugging and refactoring
    #    1   -4.508E-14  [ 1.0481009E+00  ]    -0.009300846   2.4708108E-01  [  4.3600000, 4.3600000, 4.3600000 ]
    #    2    5.462E-14  [ 1.0505473E+00  ]    -0.009752572   2.4654500E-01  [  4.3600000, 4.3600000, 4.3597350 ]
    #    3    1.095E-13  [ 1.0531431E+00  ]    -0.010230798   2.4604978E-01  [  4.3600000, 4.3600000, 4.3594551 ]
    #    4    4.419E-14  [ 1.0558985E+00  ]    -0.010737036   2.4558510E-01  [  4.3600000, 4.3600000, 4.3591594 ]
    #    5    1.186E-13  [ 1.0588245E+00  ]    -0.011272911   2.4514220E-01  [  4.3600000, 4.3600000, 4.3588468 ]
    