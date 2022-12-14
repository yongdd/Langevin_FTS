import os
import time
import numpy as np
from scipy.io import savemat
import scft

# OpenMP environment variables
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "2"  # 0, 1 or 2
os.environ["OMP_NUM_THREADS"] = "2"

# Major Simulation params
f = 0.5         # A-fraction of major BCP chain, f
eps = 1.0       # a_A/a_B, conformational asymmetry

params = {
    # "platform":"cpu-mkl",           # choose platform among [cuda, cpu-mkl]
    
    "nx":[256],         # Simulation grid numbers
    "lx":[1.38],         # Simulation box size as a_Ref * N_Ref^(1/2) unit,
                        # where "a_Ref" is reference statistical segment length
                        # and "N_Ref" is the number of segments of reference linear homopolymer chain.

    "box_is_altering":True,       # Find box size that minimizes the free energy during saddle point iteration.
    "chain_model":"continuous",   # "discrete" or "continuous" chain model
    "ds":1/200,                   # Contour step interval, which is equal to 1/N_Ref.
    "chi_n": 9.5,                 # Interaction parameter, Flory-Huggins params * N

    "segment_lengths":{         # Relative statistical segment length compared to "a_Ref.
        "A":np.sqrt(eps*eps/(eps*eps*f + (1-f))), 
        "B":np.sqrt(    1.0/(eps*eps*f + (1-f))), },

    "distinct_polymers":[{      # Distinct Polymers
        "volume_fraction":1.0,  # volume fraction of polymer chain
        "blocks":[              # AB diBlock Copolymer
            {"type":"A", "length":f, },     # A-block
            {"type":"B", "length":2*(1-f)}, # B-block
            {"type":"A", "length":f, },     # A-block
        ],},],

    "max_iter":2000,     # The maximum relaxation iterations
    "tolerance":1e-8     # Terminate iteration if the self-consistency error is less than tolerance
}

# Set initial fields
w_A = np.zeros(list(params["nx"]), dtype=np.float64)
w_B = np.zeros(list(params["nx"]), dtype=np.float64)
print("w_A and w_B are initialized to lamellar phase.")
for i in range(0,params["nx"][0]):
    w_A[i] =  np.cos(2*np.pi*i/params["nx"][0])
    w_B[i] = -np.cos(2*np.pi*i/params["nx"][0])

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
phi_A, phi_B = calculation.get_concentrations()
w_A, w_B = calculation.get_fields()

mdic = {"params":params, "dim":len(params["nx"]), "nx":params["nx"], "lx":params["lx"], "ds":params["ds"],
        "f":f, "chi_n":params["chi_n"], "epsilon":eps, "chain_model":params["chain_model"],
        "w_a":w_A, "w_b":w_B, "phi_a":phi_A, "phi_b":phi_B}
savemat("fields.mat", mdic)

# Recording first a few iteration results for debugging and refactoring
    #    1   -6.661E-16  [ 1.5458730E+00  ]    -0.004121081   2.0171039E-01  [  1.3800000 ]
    #    2    1.776E-15  [ 1.5469517E+00  ]    -0.003588250   1.5971632E-01  [  1.3803115 ]
    #    3   -1.332E-15  [ 1.5484471E+00  ]    -0.003268425   1.2739479E-01  [  1.3805301 ]
    #    4    2.442E-15  [ 1.5501717E+00  ]    -0.003080332   1.0276243E-01  [  1.3806868 ]
    #    5   -3.331E-15  [ 1.5520167E+00  ]    -0.002973330   8.4214739E-02  [  1.3808020 ]
    #    6   -4.552E-15  [ 1.5539215E+00  ]    -0.002916069   7.0449463E-02  [  1.3808893 ]
    