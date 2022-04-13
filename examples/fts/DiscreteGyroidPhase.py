# The input file produces a gyroid phase with compositional fluctuation
# For test purpose, this program stops after 200 Langevin steps
# change "langevin_max_step" to a larger number for the actual simulation
#
# -------------- Reference ------------
# M.W. Matsen, and T.M. Beardsley, Polymers 2021, 13, 2437.
# https://doi.org/10.3390/polym13152437

import sys
import os
import time
import pathlib
import numpy as np
from scipy.io import loadmat, savemat
from langevinfts import *
from find_saddle_point import *

# -------------- simulation parameters ------------
# Cuda environment variables
# os.environ["CUDA_VISIBLE_DEVICES"]= "1"
# OpenMP environment variables

os.environ["MKL_NUM_THREADS"] = "1"  # always 1
os.environ["OMP_STACKSIZE"] = "1G"
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "2"  # 0, 1 or 2

verbose_level = 1  # 1 : print at each langevin step.
                   # 2 : print at each saddle point iteration.

input_data = loadmat("GyroidInput.mat", squeeze_me=True)

# Simulation Box
nx = [64, 64, 64,]
lx = [7.31, 7.31, 7.31]

# Polymer Chain
n_contour = 90
f = 0.4
chi_n = 18.35
chain_model = "Discrete" # choose among [Gaussian, Discrete]

# Anderson Mixing
saddle_tolerance = 1e-4
saddle_max_iter = 100
am_n_comp = 1  # W+
am_max_hist= 20
am_start_error = 8e-1
am_mix_min = 0.1
am_mix_init = 0.1

# Langevin Dynamics
langevin_dt = 0.8     # langevin step interval, delta tau*N
langevin_nbar = 10000  # invariant polymerization index
langevin_max_step = 200

# -------------- initialize ------------
# choose platform among [cuda, cpu-mkl, cpu-fftw]
if "cuda" in PlatformSelector.avail_platforms():
    platform = "cuda"
else:
    platform = PlatformSelector.avail_platforms()[0]
print("platform :", platform)
factory = PlatformSelector.create_factory(platform)

# create instances
pc     = factory.create_polymer_chain(f, n_contour, chi_n, chain_model)
sb     = factory.create_simulation_box(nx, lx)
pseudo = factory.create_pseudo(sb, pc)
am     = factory.create_anderson_mixing(sb, am_n_comp,
            am_max_hist, am_start_error, am_mix_min, am_mix_init)

# standard deviation of normal noise for single segment
langevin_sigma = np.sqrt(2*langevin_dt*sb.get_n_grid()/
    (sb.get_volume()*np.sqrt(langevin_nbar)))

# random seed for MT19937
np.random.seed(5489)
# -------------- print simulation parameters ------------
print("---------- Simulation Parameters ----------")
print("Box Dimension: %d"  % (sb.get_dim()) )
print("Precision: 8")
print("chi_n: %f, f: %f, N: %d" % (pc.get_chi_n(), pc.get_f(), pc.get_n_contour()) )
print("%s chain model" % (pc.get_model_name()) )
print("Nx: %d, %d, %d" % (sb.get_nx(0), sb.get_nx(1), sb.get_nx(2)) )
print("Lx: %f, %f, %f" % (sb.get_lx(0), sb.get_lx(1), sb.get_lx(2)) )
print("dx: %f, %f, %f" % (sb.get_dx(0), sb.get_dx(1), sb.get_dx(2)) )
print("Volume: %f" % (sb.get_volume()) )

print("Invariant Polymerization Index: %d" % (langevin_nbar) )
print("Langevin Sigma: %f" % (langevin_sigma) )
print("Random Number Generator: ", np.random.RandomState().get_state()[0])

#-------------- allocate array ------------
# free end initial condition. q1 is q and q2 is qdagger.
# q1 starts from A end and q2 starts from B end.
q1_init = np.ones(sb.get_n_grid(), dtype=np.float64)
q2_init = np.ones(sb.get_n_grid(), dtype=np.float64)

print("wminus and wplus are initialized to gyroid")
w_minus = input_data["w_minus"]
w_plus = input_data["w_plus"]

# keep the level of field value
sb.zero_mean(w_plus)
sb.zero_mean(w_minus)

phi_a, phi_b = find_saddle_point(pc, sb, pseudo, am,
    q1_init, q2_init, w_plus, w_minus,
    saddle_max_iter, saddle_tolerance, verbose_level)
    
# init structure function
sf_average = np.zeros_like(np.fft.rfftn(np.reshape(w_minus, sb.get_nx()[:sb.get_dim()])),np.float64)

#------------------ run ----------------------
print("---------- Run ----------")
time_start = time.time()

print("iteration, mass error, total_partition, energy_total, error_level")
for langevin_step in range(1, langevin_max_step+1):

    print("langevin step: ", langevin_step)
    # update w_minus: predict step
    w_minus_copy = w_minus.copy()
    normal_noise = np.random.normal(0.0, langevin_sigma, sb.get_n_grid())
    lambda1 = phi_a-phi_b + 2*w_minus/pc.get_chi_n()
    w_minus += -lambda1*langevin_dt + normal_noise
    phi_a, phi_b = find_saddle_point(pc, sb, pseudo, am,
        q1_init, q2_init, w_plus, w_minus,
        saddle_max_iter, saddle_tolerance, verbose_level)

    # update w_minus: correct step
    lambda2 = phi_a-phi_b + 2*w_minus/pc.get_chi_n()
    w_minus = w_minus_copy - 0.5*(lambda1+lambda2)*langevin_dt + normal_noise
    phi_a, phi_b = find_saddle_point(pc, sb, pseudo, am,
        q1_init, q2_init, w_plus, w_minus,
        saddle_max_iter, saddle_tolerance, verbose_level)

    # calcaluate structure function
    if langevin_step % 10 == 0:
        sf_average += np.absolute(np.fft.rfftn(np.reshape(w_minus, sb.get_nx()))/sb.get_n_grid())**2

    # save structure function
    if langevin_step % 100 == 0:
        sf_average *= 10/100*sb.get_volume()*np.sqrt(langevin_nbar)/pc.get_chi_n()**2
        sf_average -= 1.0/(2*pc.get_chi_n())
        mdic = {"dim":sb.get_dim(), "nx":sb.get_nx(), "lx":sb.get_lx(),
        "N":pc.get_n_contour(), "f":pc.get_f(), "chi_n":pc.get_chi_n(),
        "chain_model":pc.get_model_name(),
        "langevin_dt":langevin_dt, "langevin_nbar":langevin_nbar,
        "structure_function":sf_average}
        savemat( "structure_function_%06d.mat" % (langevin_step), mdic)
        sf_average[:,:,:] = 0.0

    # write density and field data
    if langevin_step % 100 == 0:
        mdic = {"dim":sb.get_dim(), "nx":sb.get_nx(), "lx":sb.get_lx(),
            "N":pc.get_n_contour(), "f":pc.get_f(), "chi_n":pc.get_chi_n(),
            "chain_model":pc.get_model_name(), "n_bar":langevin_nbar,
            "random_seed":np.random.RandomState().get_state()[0],
            "w_plus":w_plus, "w_minus":w_minus, "phi_a":phi_a, "phi_b":phi_b}
        savemat( "fields_%06d.mat" % (langevin_step), mdic)

# estimate execution time
time_duration = time.time() - time_start
print("total time: %f, time per step: %f" %
    (time_duration, time_duration/langevin_max_step) )
