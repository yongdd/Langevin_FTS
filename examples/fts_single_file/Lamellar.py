# This program will produce a lamellar phase from random initial condition
# For test purpose, this program stops after 2000 Langevin steps
# change "langevin_max_step" to a larger number for the actual simulation

import os
import time
import numpy as np
from scipy.io import savemat
from langevinfts import *

def find_saddle_point(cb, pseudo, am, chi_n,
    q1_init, q2_init, w_plus, w_minus, 
    saddle_max_iter, saddle_tolerance, verbose_level):
        
    # assign large initial value for the energy and error
    energy_total = 1e20
    error_level = 1e20

    # reset Anderson mixing module
    am.reset_count()

    # saddle point iteration begins here
    for saddle_iter in range(1,saddle_max_iter+1):
        
        # for the given fields find the polymer statistics
        phi, Q = pseudo.compute_statistics(
            q1_init, q2_init,
            np.stack((w_plus+w_minus,w_plus-w_minus), axis=0))
        phi_plus = phi[0] + phi[1]
        
        # calculate output fields
        g_plus = phi_plus-1.0
        w_plus_out = w_plus + g_plus 
        cb.zero_mean(w_plus_out)

        # error_level measures the "relative distance" between the input and output fields
        old_error_level = error_level
        error_level = np.sqrt(cb.inner_product(phi_plus-1.0,phi_plus-1.0)/cb.get_volume())

        # print iteration # and error levels
        if(verbose_level == 2 or
         verbose_level == 1 and
         (error_level < saddle_tolerance or saddle_iter == saddle_max_iter)):
             
            # calculate the total energy
            energy_total  = -np.log(Q/cb.get_volume())
            energy_total += cb.inner_product(w_minus,w_minus)/chi_n/cb.get_volume()
            energy_total += chi_n/4
            energy_total -= cb.integral(w_plus)/cb.get_volume()

            # check the mass conservation
            mass_error = cb.integral(phi_plus)/cb.get_volume() - 1.0
            print("%8d %12.3E %15.7E %15.9f %15.7E" %
                (saddle_iter, mass_error, Q, energy_total, error_level))
        # conditions to end the iteration
        if error_level < saddle_tolerance:
            break
            
        # calculate new fields using simple and Anderson mixing
        am.calculate_new_fields(w_plus, w_plus_out, g_plus, old_error_level, error_level)
    return phi, Q

# -------------- simulation parameters ------------

# OpenMP environment variables
os.environ["OMP_STACKSIZE"] = "1G"
os.environ["MKL_NUM_THREADS"] = "1"  # always 1
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "2"  # 0, 1 or 2

verbose_level = 1  # 1 : print at each langevin step.
                   # 2 : print at each saddle point iteration.

# Simulation Grids and Lengths
nx = [32, 32, 32]
lx = [8.0, 8.0, 8.0]

# Polymer Chain
f = 0.5
n_segment = 16
chi_n = 20
epsilon = 1.0                # a_A/a_B, conformational asymmetry
chain_model = "Continuous"   # choose among [Continuous, Discrete]
ds = 1/n_segment             # contour step interval

# Anderson Mixing
saddle_tolerance = 1e-4
saddle_max_iter = 100
am_n_var = np.prod(nx)  # W+
am_max_hist= 20
am_start_error = 8e-1
am_mix_min = 0.1
am_mix_init = 0.1

# Langevin Dynamics
langevin_dt = 0.8        # langevin step interval, delta tau*N
langevin_nbar = 1024     # invariant polymerization index
langevin_max_step = 2000

# -------------- initialize ------------
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

# standard deviation of normal noise
langevin_sigma = np.sqrt(2*langevin_dt*cb.get_n_grid()/
    (cb.get_volume()*np.sqrt(langevin_nbar)))

## random seed for MT19937
np.random.seed(5489)
# -------------- print simulation parameters ------------
print("---------- Simulation Parameters ----------")
print("Box Dimension: %d"  % (cb.get_dim()) )
print("chi_n: %f, f: %f, N: %d" % (chi_n, f, pc.get_n_segment_total()) )
print("%s chain model" % (pc.get_model_name()) )
print("Conformational asymmetry (epsilon): %f" % (epsilon) )
print("Nx: %d, %d, %d" % (cb.get_nx(0), cb.get_nx(1), cb.get_nx(2)) )
print("Lx: %f, %f, %f" % (cb.get_lx(0), cb.get_lx(1), cb.get_lx(2)) )
print("dx: %f, %f, %f" % (cb.get_dx(0), cb.get_dx(1), cb.get_dx(2)) )
print("Volume: %f" % (cb.get_volume()) )

print("Invariant Polymerization Index: %d" % (langevin_nbar) )
print("Langevin Sigma: %f" % (langevin_sigma) )
print("Random Number Generator: ", np.random.RandomState().get_state()[0])

#-------------- allocate array ------------
# free end initial condition. q1 is q and q2 is qdagger.
# q1 starts from A end and q2 starts from B end.
q1_init = np.ones(cb.get_n_grid(), dtype=np.float64)
q2_init = np.ones(cb.get_n_grid(), dtype=np.float64)

print("w_minus and w_plus are initialized to random")
w_plus  = np.random.normal(0.0, langevin_sigma, cb.get_n_grid())
w_minus = np.random.normal(0.0, langevin_sigma, cb.get_n_grid())

# keep the level of field value
cb.zero_mean(w_plus)

phi, _ = find_saddle_point(cb, pseudo, am, chi_n,
    q1_init, q2_init, w_plus, w_minus,
    saddle_max_iter, saddle_tolerance, verbose_level)
    
# init structure function
sf_average = np.zeros_like(np.fft.rfftn(np.reshape(w_minus, cb.get_nx())),np.float64)

#------------------ run ----------------------
print("---------- Run ----------")
time_start = time.time()

print("iteration, mass error, total_partition, energy_total, error_level")
for langevin_step in range(1, langevin_max_step+1):

    print("langevin step: ", langevin_step)
    # update w_minus: predict step
    w_minus_copy = w_minus.copy()
    normal_noise = np.random.normal(0.0, langevin_sigma, cb.get_n_grid())
    lambda1 = phi[0]-phi[1] + 2*w_minus/chi_n
    w_minus += -lambda1*langevin_dt + normal_noise
    phi, _ = find_saddle_point(cb, pseudo, am, chi_n,
        q1_init, q2_init, w_plus, w_minus,
        saddle_max_iter, saddle_tolerance, verbose_level)

    # update w_minus: correct step
    lambda2 = phi[0]-phi[1] + 2*w_minus/chi_n
    w_minus = w_minus_copy - 0.5*(lambda1+lambda2)*langevin_dt + normal_noise
    phi, _ = find_saddle_point(cb, pseudo, am, chi_n,
        q1_init, q2_init, w_plus, w_minus,
        saddle_max_iter, saddle_tolerance, verbose_level)
        
    # calcaluate structure function
    if langevin_step % 10 == 0:
        sf_average += np.absolute(np.fft.rfftn(np.reshape(w_minus, cb.get_nx()))/cb.get_n_grid())**2

    # save structure function
    if langevin_step % 1000 == 0:
        sf_average *= 10/1000*cb.get_volume()*np.sqrt(langevin_nbar)/chi_n**2
        sf_average -= 1.0/(2*chi_n)
        mdic = {"dim":cb.get_dim(), "nx":cb.get_nx(), "lx":cb.get_lx(),
        "N":pc.get_n_segment_total(), "f":f, "chi_n":chi_n, "epsilon":epsilon,
        "chain_model":pc.get_model_name(),
        "dt":langevin_dt, "nbar":langevin_nbar,
        "structure_function":sf_average}
        savemat( "structure_function_%06d.mat" % (langevin_step), mdic)
        sf_average[:,:,:] = 0.0

    # write density and field data
    if langevin_step % 1000 == 0:
        mdic = {"dim":cb.get_dim(), "nx":cb.get_nx(), "lx":cb.get_lx(),
            "N":pc.get_n_segment_total(), "f":f, "chi_n":chi_n, "epsilon":epsilon,
            "chain_model":pc.get_model_name(), "nbar":langevin_nbar,
            "random_generator":np.random.RandomState().get_state()[0],
            "random_seed":np.random.RandomState().get_state()[1],
            "w_plus":w_plus, "w_minus":w_minus, "phi_a":phi[0], "phi_b":phi[1]}
        savemat( "fields_%06d.mat" % (langevin_step), mdic)

# estimate execution time
time_duration = time.time() - time_start
print("total time: %f, time per step: %f" %
    (time_duration, time_duration/langevin_max_step) )

# Recording first a few iteration results for debugging and refactoring

# w_minus and w_plus are initialized to random
#       23   -2.220E-16   5.8293236E+02     5.031220086   9.5822227E-05
# ---------- Run ----------
# iteration, mass error, total_partition, energy_total, error_level
# langevin step:  1
#       19   -3.553E-15   6.5966016E+02     5.062036019   7.6462958E-05
#        8   -2.998E-15   6.5907274E+02     5.059683923   7.3716113E-05
# langevin step:  2
#       19   -2.220E-16   7.4356322E+02     5.092392378   8.9084460E-05
#        8    4.441E-16   7.4257836E+02     5.089977265   9.0169219E-05
# langevin step:  3
#       20    2.665E-15   8.2950944E+02     5.121605195   7.2575434E-05
#        9   -3.553E-15   8.2821693E+02     5.119221963   7.0321982E-05