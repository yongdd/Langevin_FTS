import os
import time
import numpy as np
from scipy.io import savemat
from scipy.ndimage.filters import gaussian_filter
import scft

# # Major Simulation params
f = 24/90       # A-fraction of major BCP chain, f
eps = 1.0       # a_A/a_B, conformational asymmetry

params = {
    "nx":[48,48,48],            # Simulation grid numbers
    "lx":[5.74,5.74,5.74],         # Simulation box size as a_Ref * N_Ref^(1/2) unit,
                                # where "a_Ref" is reference statistical segment length
                                # and "N_Ref" is the number of segments of reference linear homopolymer chain.

    "box_is_altering":True,     # Find box size that minimizes the free energy during saddle point iteration.
    "chain_model":"discrete",   # "discrete" or "continuous" chain model
    "ds":1/90,                  # Contour step interval, which is equal to 1/N_Ref.
    "chi_n": 18.1,              # Interaction parameter, Flory-Huggins params * N

    "segment_lengths":{         # Relative statistical segment length compared to "a_Ref.
        "A":np.sqrt(eps*eps/(eps*eps*f + (1-f))), 
        "B":np.sqrt(    1.0/(eps*eps*f + (1-f))), },

    "distinct_polymers":[{      # Distinct Polymers
        "volume_fraction":1.0,  # volume fraction of polymer chain
        "blocks":[              # AB diBlock Copolymer
            {"type":"A", "length":f, }, # A-block
            {"type":"B", "length":1-f}, # B-block
        ],},],
    
    "max_iter":2000,     # The maximum relaxation iterations
    "tolerance":1e-8     # Terminate iteration if the self-consistency error is less than tolerance
}

# Set initial fields
w_A = np.zeros(list(params["nx"]), dtype=np.float64)
w_B = np.zeros(list(params["nx"]), dtype=np.float64)
print("w_A and w_B are initialized to cylindrical phase.")
n_unitcell = 3 # number of unit cell for each direction. the number of total unit cells is n_unitcell^3
sphere_positions = []
for i in range(0,n_unitcell):
    for j in range(0,n_unitcell):
        for k in range(0,n_unitcell):
            sphere_positions.append([i/n_unitcell,j/n_unitcell,k/n_unitcell])
            sphere_positions.append([(i+1/2)/n_unitcell,(j+1/2)/n_unitcell,(k+1/2)/n_unitcell])
for x,y,z in sphere_positions:
    mx, my, mz = np.round((np.array([x, y, z])*params["nx"])).astype(np.int32)
    w_A[mx,my,mz] = -1/(np.prod(params["lx"])/np.prod(params["nx"]))
w_A = gaussian_filter(w_A, sigma=np.min(params["nx"])/15, mode='wrap')

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
    #    1   -3.231E-14  [ 1.8933680E+02  ]    -0.000561236   1.3269544E+00  [  5.7400000, 5.7400000, 5.7400000 ]
    #    2   -1.125E-13  [ 1.8926963E+02  ]    -0.000298706   1.0572487E+00  [  5.7400038, 5.7400038, 5.7400038 ]
    #    3    5.329E-15  [ 1.8923283E+02  ]    -0.000164119   8.1326847E-01  [  5.7400055, 5.7400055, 5.7400055 ]
    #    4    1.910E-14  [ 1.8921242E+02  ]    -0.000095190   6.0480704E-01  [  5.7400063, 5.7400063, 5.7400063 ]
    #    5   -3.197E-14  [ 1.8920106E+02  ]    -0.000060061   4.3534872E-01  [  5.7400065, 5.7400065, 5.7400065 ]