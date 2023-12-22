#include <iostream>
#include <cmath>
#include <algorithm>

#include "Exception.h"
#include "Polymer.h"
#include "PropagatorAnalyzer.h"
#ifdef USE_CPU_MKL
#include "CpuComputationBox.h"
#include "CpuComputationContinuous.h"
#endif
#ifdef USE_CUDA
#include "CudaComputationBox.h"
#include "CudaComputationContinuous.h"
#include "CudaComputationReduceMemoryContinuous.h"
#endif

int main()
{
    try
    {
        const int II{10};
        const int M{II};
        const int NN{4};

        double q_prev[M], q_next[M];

        std::array<double,M> diff_sq;
        double error;
        double Lx, f;

        f = 0.5;
        Lx = 4.0;

        double w_a[M] = {0.822383458999126,  0.180877073118435, 0.885692320279145,
                       0.89417010799111,   0.495074990864166, 0.612975741629382,
                       0.0415795198090432, 0.353431810889399, 0.773118461366249,
                       0.474587294381635};

        double w_b[M] = {0.200002276650706, 0.592127025743285, 0.460207620078036,
                       0.435945198378862, 0.61269588805607,  0.355979618324841,
                       0.548759176402544, 0.482897565408353, 0.541501788021353,
                       0.349682106604464};

        //-------------- initialize ------------
        std::cout<< "Initializing" << std::endl;
        std::map<std::string, double> bond_lengths = {{"A",1.0}, {"B",1.0}};
        std::vector<BlockInput> blocks = 
        {
            {"A",    f, 0, 1},
            {"B",1.0-f, 1, 2},
        };

        std::vector<std::string> bc = {"absorbing", "reflecting"};

        Molecules* molecules = new Molecules("Continuous", 1.0/NN, bond_lengths);
        molecules->add_polymer(1.0, blocks, {});
        PropagatorAnalyzer* propagator_analyzer= new PropagatorAnalyzer(molecules, false);

        propagator_analyzer->display_blocks();
        propagator_analyzer->display_propagators();

        std::vector<PropagatorComputation*> solver_list;
        std::vector<std::string> solver_name_list;
        #ifdef USE_CPU_MKL
        solver_name_list.push_back("cpu-mkl");
        solver_list.push_back(new CpuComputationContinuous(new CpuComputationBox({II}, {Lx}, bc), molecules, propagator_analyzer, "realspace"));
        #endif
        #ifdef USE_CUDA
        solver_name_list.push_back("cuda");
        solver_name_list.push_back("cuda_reduce_memory_usage");
        solver_list.push_back(new CudaComputationContinuous(new CudaComputationBox({II}, {Lx}, bc), molecules, propagator_analyzer, "realspace"));
        solver_list.push_back(new CudaComputationReduceMemoryContinuous(new CudaComputationBox({II}, {Lx}, bc), molecules, propagator_analyzer, "realspace"));
        #endif

        // For each platform
        for(size_t n=0; n<solver_list.size(); n++)
        {
            PropagatorComputation* solver = solver_list[n];

            for(int i=0; i<M; i++)
                q_next[i] = 0.0;

            //---------------- run --------------------
            std::cout<< std::endl << "Running Pseudo: " << solver_name_list[n] << std::endl;
            solver->compute_statistics({{"A",w_a},{"B",w_b}},{});

            Polymer& pc = molecules->get_polymer(0);
            solver->get_chain_propagator(q_next, 0, 1, 2, pc.get_block(1,2).n_segment);
            if (n > 0)
            {
                //--------------- check --------------------
                std::cout<< "Checking"<< std::endl;
                std::cout<< "If error is less than 1.0e-7, it is ok!" << std::endl;
                std::cout<< "i: q (" << solver_name_list[n-1] << "), (" << solver_name_list[n] << ")" << std::endl;
                for(int i=0; i<M; i++)
                    std::cout << i << ": " << q_prev[i] << ", " << q_next[i] << std::endl;

                for(int i=0; i<M; i++)
                    diff_sq[i] = pow(q_prev[i] - q_next[i],2);
                error = sqrt(*std::max_element(diff_sq.begin(), diff_sq.end()));

                std::cout<< "Propagator error: "<< error << std::endl;
                if (!std::isfinite(error) || error > 1e-7)
                    return -1;
            }

            for(int i=0; i<M; i++)
                q_prev[i] = q_next[i];
        }
        return 0;
    }
    catch(std::exception& exc)
    {
        std::cout << exc.what() << std::endl;
        return -1;
    }
}
