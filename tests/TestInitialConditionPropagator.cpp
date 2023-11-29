#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <array>
#include <chrono>

#include "Exception.h"
#include "ComputationBox.h"
#include "Polymer.h"
#include "Molecules.h"
#include "Pseudo.h"
#include "AndersonMixing.h"
#include "AbstractFactory.h"
#include "PlatformSelector.h"

int main()
{
    try
    {
        // string to output file and print stream
        std::streamsize default_precision = std::cout.precision();

        // -------------- initialize ------------
        // platform type, [cuda, cpu-mkl]

        std::vector<int> nx = {16,16,16};
        std::vector<double> lx = {4.0, 4.0, 4.0};
        const int N = 100;
        double ds = 1.0/N;

        std::map<std::string, double> bond_lengths = {{"A",1.0}};
        std::vector<BlockInput> block_inputs = {{"A", 1.0, 0, 1}};

        // from a vertex index to a grafting point
        // following map means that vertex 
        std::map<int, std::string> chain_end_to_q_init = {{0,"G"}};
        const int M = nx[0]*nx[1]*nx[2];

        bool reduce_memory_usage=false;

        //-------------- allocate array ------------
        double w[M];
        double q_init[M];
        double q_out[M];

        // choose platform
        std::vector<std::string> chain_models = {"Continuous", "Discrete"};
        std::vector<std::string> avail_platforms = PlatformSelector::avail_platforms();
        std::vector<bool> aggregate_propagator_computations = {false, true};
        for(std::string chain_model : chain_models)
        {
            std::vector<double> x_square_total_list;
            for(std::string platform : avail_platforms)
            {
                for(bool aggregate_propagator_computation : aggregate_propagator_computations)
                {
                    AbstractFactory *factory = PlatformSelector::create_factory(platform, reduce_memory_usage);
                    // factory->display_info();

                    // create instances and assign to the variables of base classes for the dynamic binding
                    ComputationBox *cb = factory->create_computation_box(nx, lx);
                    Molecules* molecules        = factory->create_molecule_information(chain_model, ds, bond_lengths, aggregate_propagator_computation);
                    molecules->add_polymer(1.0, block_inputs, chain_end_to_q_init);
                    Pseudo *pseudo     = factory->create_pseudo(cb, molecules);

                    // -------------- print simulation parameters ------------
                    std::cout << std::setprecision(default_precision);
                    std::cout << "Chain Model: " << molecules->get_model_name() << std::endl;
                    std::cout << "Platform: " << platform << std::endl;
                    std::cout << "Using Aggregation: " << aggregate_propagator_computation << std::endl;

                    // // display branches
                    // molecules->display_blocks();
                    // molecules->display_sub_propagators();

                    for(int i=0; i<M; i++)
                    {
                        w[i] = 0.0;
                        q_init[i] = 0.0;
                    }
                    q_init[0] = 1.0/cb->get_dv(0);

                    // keep the level of field value
                    cb->zero_mean(w);

                    // for the given fields find the polymer statistics
                    pseudo->compute_statistics({{"A",w}},{{"G", q_init}});

                    for(int n=20; n<=N; n+=20)
                    {   
                                                   //output, p, v ,u, n
                        pseudo->get_chain_propagator(q_out, 0, 0, 1, n);

                        double sum = 0.0;
                        double x_square = 0.0;
                        double xx, yy, zz;
                        int idx;
                        for(int i=0; i<nx[0]; i++)
                        {
                            xx = cb->get_dx(0)*std::min(i, nx[0]-i);
                            for(int j=0; j<nx[1]; j++)
                            {
                                yy = cb->get_dx(1)*std::min(j, nx[1]-j);
                                for(int k=0; k<nx[2]; k++)
                                {
                                    zz = cb->get_dx(2)*std::min(k, nx[2]-k);
                                    idx = i*nx[1]*nx[2] + j*nx[2] + k;
                                    x_square += q_out[idx]*cb->get_dv(idx)*(xx*xx + yy*yy + zz*zz);
                                    sum += q_out[idx]*cb->get_dv(idx);
                                }
                            }
                        }
                        if (chain_model == "Continuous")
                        {
                            x_square *= ((double) N)/((double) n);
                            std::cout << "n, <x^2>N/n: " << n << ", " << x_square << std::endl;
                        }
                        else if (chain_model == "Discrete")
                        {
                            x_square *= ((double) N)/((double) n-1);
                            std::cout << "n, <x^2>N/(n-1): " << n << ", " << x_square << std::endl;
                        }
                        if ( std::abs(x_square-1.0) > 1e-2)
                        {
                            std::cout << "The error of x_square (" << x_square-1.0 << ") is not small." << std::endl;
                            return -1;
                        }
                    }
                    
                    delete factory;
                    delete cb;
                    delete molecules;
                    delete pseudo;
                }
            }
        }
        return 0;
    }
    catch(std::exception& exc)
    {
        std::cout << exc.what() << std::endl;
        return -1;
    }
}