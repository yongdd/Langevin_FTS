#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <array>
#include <chrono>

#include "Exception.h"
#include "ParamParser.h"
#include "BranchedPolymerChain.h"
#include "ComputationBox.h"
#include "PseudoBranched.h"
#include "AndersonMixing.h"
#include "AbstractFactory.h"
#include "PlatformSelector.h"

int main()
{
    try
    {
        // math constants
        const double PI = 3.14159265358979323846;
        // chrono timer
        std::chrono::system_clock::time_point chrono_start, chrono_end;
        std::chrono::duration<double> time_duration;

        // QQ = total partition function
        double QQ_a, QQ_b, energy_total;
        // error_level = variable to check convergence of the iteration
        double error_level, old_error_level;
        // input and output fields, xi is temporary storage for pressures
        double *w, *w_out, *w_diff;  // n_comp * MM
        double *xi, *w_plus, *w_minus; // MM
        // initial value of q, q_dagger
        double *q1_init, *q2_init;
        // segment concentration
        double *phi, *phi_a, *phi_b, *phitot;

        // string to output file and print stream
        std::ofstream print_stream;
        std::stringstream ss;
        std::string print_file_name;
        // temp
        int idx;
        double sum;

        // -------------- initialize ------------
        // platform type, [cuda, cpu-mkl]
        
        int max_scft_iter = 20;
        double tolerance = 1e-9;

        double frac_a = 0.5;

        double chi_n = 3.0;
        std::vector<int> nx = {263};
        std::vector<double> lx = {5.0};
        std::string chain_model = "Discrete";  // choose among [Continuous, Discrete]
        double ds = 1.0/50;

        int am_n_var = 2*nx[0];  // A and B
        int am_max_hist= 20;
        double am_start_error = 8e-1;
        double am_mix_min = 0.1;
        double am_mix_init = 0.1;

        // choose platform
        std::vector<std::string> avail_platforms = PlatformSelector::avail_platforms();
        for(std::string platform : avail_platforms){
            AbstractFactory *factory = PlatformSelector::create_factory(platform,chain_model);
            factory->display_info();

            // create instances and assign to the variables of base classes for the dynamic binding
            ComputationBox *cb  = factory->create_computation_box(nx, lx);
            BranchedPolymerChain *pc_a = factory->create_polymer_chain(ds, {{"A",1.0}}, {"A"}, {1.0}, {0}, {1}, {});
            BranchedPolymerChain *pc_b = factory->create_polymer_chain(ds, {{"B",1.0}}, {"B"}, {1.0}, {0}, {1}, {});
            PseudoBranched *pseudo_a   = factory->create_pseudo(cb, pc_a);
            PseudoBranched *pseudo_b   = factory->create_pseudo(cb, pc_b);
            AndersonMixing *am = factory->create_anderson_mixing(am_n_var,
                                am_max_hist, am_start_error, am_mix_min, am_mix_init);

            // -------------- print simulation parameters ------------
            std::cout<< "---------- Simulation Parameters ----------" << std::endl;
            std::cout << "Box Dimension: " << cb->get_dim() << std::endl;
            std::cout << "chi_n, frac_a, N_a, N_b: " << chi_n << " " << frac_a << " " << pc_a->get_n_segment_total() << " " << pc_b->get_n_segment_total() << std::endl;
            std::cout << "Nx: " << cb->get_nx(0) << " " << cb->get_nx(1) << " " << cb->get_nx(2) << std::endl;
            std::cout << "Lx: " << cb->get_lx(0) << " " << cb->get_lx(1) << " " << cb->get_lx(2) << std::endl;
            std::cout << "dx: " << cb->get_dx(0) << " " << cb->get_dx(1) << " " << cb->get_dx(2) << std::endl;
            sum = 0.0;
            for(int i=0; i<cb->get_n_grid(); i++)
                sum += cb->get_dv(i);
            std::cout << "volume, sum(dv):  " << cb->get_volume() << " " << sum << std::endl;

            //-------------- allocate array ------------
            w       = new double[cb->get_n_grid()*2];
            w_out   = new double[cb->get_n_grid()*2];
            w_diff  = new double[cb->get_n_grid()*2];
            xi      = new double[cb->get_n_grid()];
            phi     = new double[cb->get_n_grid()*2];
            phitot  = new double[cb->get_n_grid()];
            w_plus  = new double[cb->get_n_grid()];
            w_minus = new double[cb->get_n_grid()];
            q1_init = new double[cb->get_n_grid()];
            q2_init = new double[cb->get_n_grid()];

            phi_a = &phi[0];
            phi_b = &phi[cb->get_n_grid()];
            //-------------- setup fields ------------
            //call random_number(phi_a)
            //   phi_a = reshape( phi_a, (/ x_hi-x_lo+1,y_hi-y_lo+1,z_hi-z_lo+1 /), order = (/ 3, 2, 1 /))
            //   call random_number(phi_a(:,:,z_lo))
            //   do k=z_lo,z_hi
            //     phi_a(:,:,k) = phi_a(:,:,z_lo)
            //   end do

            std::cout<< "w_a and w_b are initialized to a given test fields." << std::endl;
            for(int i=0; i<cb->get_nx(0); i++)
                for(int j=0; j<cb->get_nx(1); j++)
                    for(int k=0; k<cb->get_nx(2); k++)
                    {
                        idx = i*cb->get_nx(1)*cb->get_nx(2) + j*cb->get_nx(2) + k;
                        phi_a[idx]= cos(2*(cb->get_nx(2)-k)/cb->get_nx(2)*PI);
                    }

            for(int i=0; i<cb->get_n_grid(); i++)
            {
                phi_b[i] = 1.0 - phi_a[i];
                w[i]                  = chi_n*phi_b[i];
                w[i+cb->get_n_grid()] = chi_n*phi_a[i];
            }

            // keep the level of field value
            cb->zero_mean(&w[0]);
            cb->zero_mean(&w[cb->get_n_grid()]);

            // free end initial condition. q1 is q and q2 is qdagger.
            // q1 starts from A end and q2 starts from B end.
            for(int i=0; i<cb->get_n_grid(); i++)
            {
                q1_init[i] = 1.0;
                q2_init[i] = 1.0;
            }
            //print_stream->close();

            // assign large initial value for the energy and error
            energy_total = 1.0e20;
            error_level = 1.0e20;

            //------------------ run ----------------------
            std::cout<< "---------- Run ----------" << std::endl;
            std::cout<< "iteration, mass error, total_partition_a, _b, energy_total, error_level" << std::endl;
            chrono_start = std::chrono::system_clock::now();
            // iteration begins here
            for(int iter=0; iter<max_scft_iter; iter++)
            {
                // for the given fields find the polymer statistics
                pseudo_a->compute_statistics({}, {{"A",&w[0]}},                phi_a, QQ_a);
                pseudo_b->compute_statistics({}, {{"B",&w[cb->get_n_grid()]}}, phi_b, QQ_b);

                for(int i=0; i<cb->get_n_grid(); i++)
                {
                    phi_a[i] = phi_a[i]*frac_a;
                    phi_b[i] = phi_b[i]*(1.0-frac_a);
                }

                // calculate the total energy
                for(int i=0; i<cb->get_n_grid(); i++)
                {
                    w_minus[i] = (w[i]-w[i + cb->get_n_grid()])/2;
                    w_plus[i]  = (w[i]+w[i + cb->get_n_grid()])/2;
                }
                energy_total  = -frac_a*log(QQ_a/cb->get_volume());
                energy_total -= (1.0-frac_a)*log(QQ_b/cb->get_volume());
                energy_total += cb->inner_product(w_minus,w_minus)/chi_n/cb->get_volume();
                energy_total += cb->integral(w_plus)/cb->get_volume();
                //energy_total += cb->inner_product(ext_w_minus,ext_w_minus)/chi_b/cb->get_volume();

                for(int i=0; i<cb->get_n_grid(); i++)
                {
                    // calculate pressure field for the new field calculation, the method is modified from Fredrickson's
                    xi[i] = 0.5*(w[i]+w[i+cb->get_n_grid()]-chi_n);
                    // calculate output fields
                    w_out[i]                  = chi_n*phi_b[i] + xi[i];
                    w_out[i+cb->get_n_grid()] = chi_n*phi_a[i] + xi[i];
                }
                cb->zero_mean(&w_out[0]);
                cb->zero_mean(&w_out[cb->get_n_grid()]);

                // error_level measures the "relative distance" between the input and output fields
                old_error_level = error_level;
                for(int i=0; i<2*cb->get_n_grid(); i++)
                    w_diff[i] = w_out[i]- w[i];
                error_level = sqrt(cb->multi_inner_product(2,w_diff,w_diff)/
                                (cb->multi_inner_product(2,w,w)+1.0));

                // print iteration # and error levels and check the mass conservation
                sum = (cb->integral(phi_a) + cb->integral(phi_b))/cb->get_volume() - 1.0;
                std::cout<< std::setw(8) << iter;
                std::cout<< std::setw(13) << std::setprecision(3) << std::scientific << sum ;
                std::cout<< std::setw(17) << std::setprecision(7) << std::scientific << QQ_a;
                std::cout<< std::setw(17) << std::setprecision(7) << std::scientific << QQ_b;
                std::cout<< std::setw(15) << std::setprecision(9) << std::fixed << energy_total;
                std::cout<< std::setw(15) << std::setprecision(9) << std::fixed << error_level << std::endl;

                // conditions to end the iteration
                if(error_level < tolerance) break;
                // calculte new fields using simple and Anderson mixing
                am->calculate_new_fields(w, w_out, w_diff, old_error_level, error_level);
            }

            // estimate execution time
            chrono_end = std::chrono::system_clock::now();
            time_duration = chrono_end - chrono_start;
            std::cout<< "total time, time per step: ";
            std::cout<< time_duration.count() << " " << time_duration.count()/max_scft_iter << std::endl;

            //------------- finalize -------------
            delete[] w;
            delete[] w_out;
            delete[] w_diff;
            delete[] xi;
            delete[] phi;
            delete[] phitot;
            delete[] w_plus;
            delete[] w_minus;
            delete[] q1_init;
            delete[] q2_init;

            delete pc_a;
            delete pc_b;
            delete cb;
            delete pseudo_a;
            delete pseudo_b;
            delete am;
            delete factory;
            if (std::isnan(error_level) || std::abs(error_level-0.000823320) > 1e-7)
                return -1;
        }
        return 0;
    }
    catch(std::exception& exc)
    {
        std::cout << exc.what() << std::endl;
        return -1;
    }
}
