import os
import numpy as np
import time
from scipy.io import loadmat, savemat
from scipy.ndimage.filters import gaussian_filter
from langevinfts import *
from find_saddle_point import *

# -------------- initialize ------------

# OpenMP environment variables
os.environ["MKL_NUM_THREADS"] = "1"  # always 1
os.environ["OMP_STACKSIZE"] = "1G"
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "2"  # 0, 1 or 2
os.environ["OMP_NUM_THREADS"] = "2"

max_scft_iter = 2000
tolerance = 1e-8

# Major Simulation Parameters
f = 0.25          # A-fraction, f
n_segment = 100   # segment number, N
chi_n = 25        # Flory-Huggins Parameters * N
epsilon = 2.0     # a_A/a_B, conformational asymmetry, default = 1.0
nx = [64,64,32]   # grid numbers
lx = [7.0,7.0,4.0]   # as aN^(1/2) unit, a = sqrt(f*a_A^2 + (1-f)*a_B^2)
ds = 1/n_segment      # contour step interval
chain_model = "Continuous"    # choose among [Continuous, Discrete]

# Anderson mixing
am_n_var = 2*np.prod(nx)+len(lx)  # w_a (w[0]) and w_b (w[1]) + lx
am_max_hist= 20                   # maximum number of history
am_start_error = 1e-2             # when switch to AM from simple mixing
am_mix_min = 0.1                  # minimum mixing rate of simple mixing
am_mix_init = 0.1                 # initial mixing rate of simple mixing

# calculate chain parameters
# a : statistical segment length, N: n_segment
# a_sq_n = [a_A^2 * N, a_B^2 * N]
a_sq_n = [epsilon*epsilon/(f*epsilon*epsilon + (1.0-f)),
            1.0/(f*epsilon*epsilon + (1.0-f))]
N_pc = [int(f*n_segment),int((1-f)*n_segment)]

# choose platform among [cuda, cpu-mkl]
if "cuda" in PlatformSelector.avail_platforms():
    platform = "cuda"
else:
    platform = PlatformSelector.avail_platforms()[0]
print("platform :", platform)
factory = PlatformSelector.create_factory(platform, chain_model)

# create instances
pc     = factory.create_polymer_chain(N_pc, np.sqrt(a_sq_n), ds)
cb     = factory.create_computation_box(nx, lx)
pseudo = factory.create_pseudo(cb, pc)
am     = factory.create_anderson_mixing(am_n_var,
            am_max_hist, am_start_error, am_mix_min, am_mix_init)

# -------------- print simulation parameters ------------
print("---------- Simulation Parameters ----------")
print("Box Dimension: %d" % (cb.get_dim()))
print("chi_n: %f, f: %f, N: %d" % (chi_n, f, pc.get_n_segment_total()) )
print("%s chain model" % (pc.get_model_name()) )
print("Conformational asymmetry (epsilon): %f" % (epsilon) )
print("Nx: %d, %d, %d" % (cb.get_nx(0), cb.get_nx(1), cb.get_nx(2)) )
print("Lx: %f, %f, %f" % (cb.get_lx(0), cb.get_lx(1), cb.get_lx(2)) )
print("dx: %f, %f, %f" % (cb.get_dx(0), cb.get_dx(1), cb.get_dx(2)) )
print("Volume: %f" % (cb.get_volume()) )

#-------------- allocate array ------------
# free end initial condition. q1 is q and q2 is qdagger.
# q1 starts from A end and q2 starts from B end.
w       = np.zeros([2]+list(cb.get_nx()), dtype=np.float64)
q1_init = np.ones (    cb.get_n_grid(),   dtype=np.float64)
q2_init = np.ones (    cb.get_n_grid(),   dtype=np.float64)

# Initial Fields
print("w_A and w_B are initialized to Sigma phase.")
# [Ref: https://doi.org/10.3390/app2030654]
sphere_positions = [[0.00,0.00,0.00],[0.50,0.50,0.50], #A
[0.40,0.40,0.00],[0.60,0.60,0.00],[0.10,0.90,0.50],[0.90,0.10,0.50], #B
[0.13,0.46,0.00],[0.46,0.13,0.00],[0.54,0.87,0.00],[0.87,0.54,0.00], #C
[0.04,0.63,0.50],[0.63,0.04,0.50],[0.37,0.96,0.50],[0.96,0.37,0.50], #C
[0.07,0.74,0.00],[0.74,0.07,0.00],[0.26,0.93,0.00],[0.93,0.26,0.00], #D
[0.24,0.43,0.50],[0.43,0.24,0.50],[0.57,0.76,0.50],[0.77,0.56,0.50], #D
[0.18,0.18,0.25],[0.82,0.82,0.25],[0.32,0.68,0.25],[0.68,0.32,0.25], #E
[0.18,0.18,0.75],[0.82,0.82,0.75],[0.32,0.68,0.75],[0.68,0.32,0.75]] #E

for x,y,z in sphere_positions:
    mx, my, mz = np.round((np.array([x, y, z])*cb.get_nx())).astype(np.int32)
    w[0,mx,my,mz] = -1/np.prod(cb.get_dx())
w[0] = gaussian_filter(w[0], sigma=np.min(cb.get_nx())/15, mode='wrap')
w = np.reshape(w, [2, cb.get_n_grid()])

# keep the level of field value
cb.zero_mean(w[0])
cb.zero_mean(w[1])

#------------------ run ----------------------
print("---------- Run ----------")
time_start = time.time()

phi, Q, energy_total = find_saddle_point(pc, cb, pseudo, am, lx, chi_n,
    q1_init, q2_init, w, max_scft_iter, tolerance, is_box_altering=True)

# estimate execution time
time_duration = time.time() - time_start
print("total time: %f " % time_duration)

# save final results
mdic = {"dim":cb.get_dim(), "nx":cb.get_nx(), "lx":cb.get_lx(),
        "N":pc.get_n_segment_total(), "f":f, "chi_n":chi_n, "epsilon":epsilon,
        "chain_model":chain_model, "w_a":w[0], "w_b":w[1], "phi_a":phi[0], "phi_b":phi[1]}
savemat("fields.mat", mdic)