import numpy as np
from langevinfts import *

def find_saddle_point(pc, sb, pseudo, am, chi_n,
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
        phi, Q = pseudo.find_phi(
            q1_init, q2_init,
            np.stack((w_plus+w_minus,w_plus-w_minus), axis=0))
        phi = phi.reshape(pc.get_n_block(),sb.get_n_grid())
        phi_plus = phi[0] + phi[1]
        
        # calculate output fields
        g_plus = phi_plus-1.0
        w_plus_out = w_plus + g_plus 
        sb.zero_mean(w_plus_out)

        # error_level measures the "relative distance" between the input and output fields
        old_error_level = error_level
        error_level = np.sqrt(sb.inner_product(phi_plus-1.0,phi_plus-1.0)/sb.get_volume())

        # print iteration # and error levels
        if(verbose_level == 2 or
         verbose_level == 1 and
         (error_level < saddle_tolerance or saddle_iter == saddle_max_iter)):
             
            # calculate the total energy
            energy_old = energy_total
            energy_total  = -np.log(Q/sb.get_volume())
            energy_total += sb.inner_product(w_minus,w_minus)/chi_n/sb.get_volume()
            energy_total += chi_n/4
            energy_total -= sb.integral(w_plus)/sb.get_volume()

            # check the mass conservation
            mass_error = sb.integral(phi_plus)/sb.get_volume() - 1.0
            print("%8d %12.3E %15.7E %15.9f %15.7E" %
                (saddle_iter, mass_error, Q, energy_total, error_level))
        # conditions to end the iteration
        if error_level < saddle_tolerance:
            break
            
        # calculte new fields using simple and Anderson mixing
        am.caculate_new_fields(w_plus, w_plus_out, g_plus, old_error_level, error_level)
    return phi, Q
