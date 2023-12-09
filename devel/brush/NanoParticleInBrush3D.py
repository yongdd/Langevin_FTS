import os
import time
import numpy as np
from scipy.io import savemat
# import scft_brush as scft
import scft

# OpenMP environment variables
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"  # 0, 1
os.environ["OMP_NUM_THREADS"] = "2"  # 1 ~ 4

params = {
    # "platform":"cuda",           # choose platform among [cuda, cpu-mkl]
    
    "nx":[240,200,200],          # Simulation grid numbers
    "lx":[12,10,10],           # Simulation box size as a_Ref * N_Ref^(1/2) unit,
                              # where "a_Ref" is reference statistical segment length
                              # and "N_Ref" is the number of segments of reference linear homopolymer chain.
                              
    "reduce_gpu_memory_usage":False, # Reduce gpu memory usage by storing propagators in main memory instead of gpu memory.
    "box_is_altering":False,    # Find box size that minimizes the free energy during saddle point iteration.
    "chain_model":"continuous",   # "discrete" or "continuous" chain model
    "ds":1/30,                 # Contour step interval, which is equal to 1/N_Ref.

    "segment_lengths":{         # Relative statistical segment length compared to "a_Ref.
        "A":1.0},

    "distinct_polymers":[  # Distinct Polymers
        {   # A Grafted Brush
            "volume_fraction":1.0,
            "blocks":[{"type":"A", "length":1.0, "v":0, "u":1}],
            "initial_conditions":{0:"G"}},
        ],

    "optimizer":{
        "name":"am",            # Anderson Mixing
        "max_hist":20,          # Maximum number of history
        "start_error":1e-2,     # When switch to AM from simple mixing
        "mix_min":0.1,         # Minimum mixing rate of simple mixing
        "mix_init":0.1,        # Initial mixing rate of simple mixing
    },

    "max_iter":2000,     # The maximum relaxation iterations
    "tolerance":1e-8     # Terminate iteration if the self-consistency error is less than tolerance
}

# Target density profile for film
T = 1.0
t = 0.4
T_mask = 1.0
L = params["lx"][0] - 2*T_mask
dx = params["lx"][0]/params["nx"][0]
I_range = round(L/dx)-3
offset = round(T_mask/dx)+1
offset_grafting = np.max([round((T_mask+0.05)/dx), round((T_mask)/dx)+1])

phi_target = np.zeros(params["nx"])
for i in range(I_range):
    phi_target[i+offset,:,:] = 1.0
    phi_target[params["nx"][0]-offset-i-1,:,:] = 1.0
print(phi_target[:30,0,0])
print(phi_target[-30:,0,0])

# Set initial fields
w_A = np.zeros(list(params["nx"]), dtype=np.float64)
print("w_A and w_B are initialized to lamellar phase.")
for i in range(I_range):
    w_A[i+offset,:,:] = -np.cos(2*np.pi*i/I_range)
# w_A = np.random.normal(0.0, 0.1, params["nx"])

# Initial condition of q (grafting point)
q_init = {"G":np.zeros(list(params["nx"]), dtype=np.float64)}
q_init["G"][offset_grafting,:,:] = 1.0/dx
q_init["G"][params["nx"][0]-offset_grafting-1,:,:] = 1.0/dx

# Mask for Nano Particle
q_mask = np.ones(params["nx"])
nano_particle_radius = 0.7

x = np.linspace(0.0-T-1.0, params["lx"][0]-T-1.0,      num=params["nx"][0], endpoint=False)
y = np.linspace(-params["lx"][1]/2, params["lx"][1]/2, num=params["nx"][1], endpoint=False)
z = np.linspace(-params["lx"][2]/2, params["lx"][2]/2, num=params["nx"][2], endpoint=False)
xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')

q_mask[np.sqrt(xv**2+yv**2+zv**2) < nano_particle_radius] = 0.0
q_mask[np.isclose(phi_target, 0.0)] = 0.0
# print(q_mask[:,0])
# print(q_init["G"][:])
print(q_mask[20:40,0])
print(q_init["G"][20:40,0])

q_mask = q_mask*np.flip(q_mask, axis=0)
phi_target = phi_target*q_mask

# Initialize calculation
calculation = scft.SCFT(params=params)

# Set a timer
time_start = time.time()

# Run
calculation.run(initial_fields={"A": w_A}, q_init=q_init, q_mask=q_mask)

# Estimate execution time
time_duration = time.time() - time_start
print("total time: %f " % time_duration)

# Save final results
calculation.save_results("fields.mat", q_mask=q_mask)
    