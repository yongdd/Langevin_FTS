#include <complex>
#include <thrust/reduce.h>
#include <iostream>
#include "CudaComputationBox.h"
#include "CudaSolverReduceMemoryDiscrete.h"

CudaSolverReduceMemoryDiscrete::CudaSolverReduceMemoryDiscrete(
    ComputationBox *cb,
    Molecules *molecules,
    PropagatorsAnalyzer *propagators_analyzer)
    : Solver(cb, molecules, propagators_analyzer)
{
    try
    {
        const int M = cb->get_n_grid();
        const int N_GPUS = CudaCommon::get_instance().get_n_gpus();
        this->propagators_analyzer = propagators_analyzer;

        // Create streams
        for(int gpu=0; gpu<N_GPUS; gpu++)
        {
            gpu_error_check(cudaSetDevice(gpu));
            gpu_error_check(cudaStreamCreate(&streams[gpu][0])); // for kernel execution
            gpu_error_check(cudaStreamCreate(&streams[gpu][1])); // for memcpy
        }
        this->propagator_solver = new CudaPseudo(cb, molecules, streams, true);

        // Allocate memory for propagators
        if( propagators_analyzer->get_essential_propagator_codes().size() == 0)
            throw_with_line_number("There is no propagator code. Add polymers first.");
        for(const auto& item: propagators_analyzer->get_essential_propagator_codes())
        {
             // There are N segments

             // Example (N==5)
             // O--O--O--O--O
             // 0  1  2  3  4

             // Legend)
             // -- : full bond
             // O  : full segment

            std::string key = item.first;
            int max_n_segment = item.second.max_n_segment;

            propagator_size[key] = max_n_segment;
            propagator[key] = new double*[max_n_segment];
            // Allocate pinned memory for device overlapping
            for(int i=0; i<propagator_size[key]; i++)
                gpu_error_check(cudaMallocHost((void**)&propagator[key][i], sizeof(double)*M));

            #ifndef NDEBUG
            propagator_finished[key] = new bool[max_n_segment];
            for(int i=0; i<max_n_segment;i++)
                propagator_finished[key][i] = false;
            #endif
        }

        // Allocate memory for propagator_junction, which contain propagator at junction of discrete chain
        for(const auto& item: propagators_analyzer->get_essential_propagator_codes())
        {
            std::string key = item.first;
            propagator_junction[key] = nullptr;
            gpu_error_check(cudaMallocHost((void**)&propagator_junction[key], sizeof(double)*M));
        }

        // Allocate memory for concentrations
        if( propagators_analyzer->get_essential_blocks().size() == 0)
            throw_with_line_number("There is no block. Add polymers first.");
        for(const auto& item: propagators_analyzer->get_essential_blocks())
        {
            phi_block[item.first] = nullptr;
            gpu_error_check(cudaMallocHost((void**)&phi_block[item.first], sizeof(double)*M));
        }

        // Total partition functions for each polymer
        single_polymer_partitions = new double[molecules->get_n_polymer_types()];

        // Remember one segment for each polymer chain to compute total partition function
        int current_p = 0;
        for(const auto& block: phi_block)
        {
            const auto& key = block.first;
            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            // Skip if already found one segment
            if (p != current_p)
                continue;

            int n_aggregated;
            int n_segment_offset    = propagators_analyzer->get_essential_block(key).n_segment_offset;
            int n_segment_original  = propagators_analyzer->get_essential_block(key).n_segment_original;
            std::string monomer_type = propagators_analyzer->get_essential_block(key).monomer_type;

            // Contains no '['
            if (dep_u.find('[') == std::string::npos)
                n_aggregated = 1;
            else
                n_aggregated = propagators_analyzer->get_essential_block(key).v_u.size();

            single_partition_segment.push_back(std::make_tuple(
                p,
                propagator[dep_v][n_segment_original-n_segment_offset-1],  // q
                propagator[dep_u][0],                                      // q_dagger
                monomer_type,       
                n_aggregated                   // how many propagators are aggregated
                ));
            current_p++;
        }

       // Find propagators and bond length for each segment to prepare stress computation
        for(const auto& block: phi_block)
        {
            const auto& key = block.first;
            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            const int N           = propagators_analyzer->get_essential_block(key).n_segment_allocated;
            const int N_OFFSET    = propagators_analyzer->get_essential_block(key).n_segment_offset;
            const int N_ORIGINAL  = propagators_analyzer->get_essential_block(key).n_segment_original;

            double **q_1 = propagator[dep_v];    // dependency v
            double **q_2 = propagator[dep_u];    // dependency u

            auto& _block_stress_info_key = block_stress_info[key];

            // Find propagators and bond length
            for(int n=0; n<=N; n++)
            {
                double *propagator_v = nullptr;
                double *propagator_u = nullptr;
                bool is_half_bond_length = false;

                // At v
                if (n + N_OFFSET == N_ORIGINAL)
                {
                    if (propagators_analyzer->get_essential_propagator_code(dep_v).deps.size() == 0) // if v is leaf node, skip
                    {
                        _block_stress_info_key.push_back(std::make_tuple(propagator_v, propagator_u, is_half_bond_length));
                        continue;
                    }
                    
                    propagator_v = propagator_junction[dep_v];
                    propagator_u = q_2[N-1];
                    is_half_bond_length = true;
                }
                // At u
                else if (n + N_OFFSET == 0){
                    if (propagators_analyzer->get_essential_propagator_code(dep_u).deps.size() == 0) // if u is leaf node, skip
                    {
                        _block_stress_info_key.push_back(std::make_tuple(propagator_v, propagator_u, is_half_bond_length));
                        continue;
                    }

                    propagator_v = q_1[N_ORIGINAL-1];
                    propagator_u = propagator_junction[dep_u];
                    is_half_bond_length = true;
                }
                // At aggregation junction
                else if (n == 0)
                {
                    _block_stress_info_key.push_back(std::make_tuple(propagator_v, propagator_u, is_half_bond_length));
                    continue;
                }
                // Within the blocks
                else
                {
                    propagator_v = q_1[N_ORIGINAL-N_OFFSET-n-1];
                    propagator_u = q_2[n-1];
                    is_half_bond_length = false;
                }
                _block_stress_info_key.push_back(std::make_tuple(propagator_v, propagator_u, is_half_bond_length));
            }
        }

        // Total partition functions for each solvent
        single_solvent_partitions = new double[molecules->get_n_solvent_types()];

        // Concentrations for each solvent
        for(int s=0;s<molecules->get_n_solvent_types();s++)
            phi_solvent.push_back(new double[M]);

        // Create scheduler for computation of propagator
        sc = new Scheduler(propagators_analyzer->get_essential_propagator_codes(), N_SCHEDULER_STREAMS); 

        gpu_error_check(cudaSetDevice(0));
        // Allocate memory for pseudo-spectral: advance_propagator()
        gpu_error_check(cudaMalloc((void**)&d_propagator_sub_dep[0], sizeof(double)*M)); // for prev
        gpu_error_check(cudaMalloc((void**)&d_propagator_sub_dep[1], sizeof(double)*M)); // for next

        for(int gpu=0; gpu<N_GPUS; gpu++)
        {
            gpu_error_check(cudaSetDevice(gpu));
            // Allocate memory for propagator computation
            gpu_error_check(cudaMalloc((void**)&d_q_one[gpu][0], sizeof(double)*M)); // for prev
            gpu_error_check(cudaMalloc((void**)&d_q_one[gpu][1], sizeof(double)*M)); // for next
            // Allocate memory for pseudo-spectral: advance_propagator()
            gpu_error_check(cudaMalloc((void**)&d_q_mask[gpu], sizeof(double)*M));
        }

        gpu_error_check(cudaSetDevice(0));
        gpu_error_check(cudaMalloc((void**)&d_q_two[0], sizeof(double)*2*M)); // for prev
        gpu_error_check(cudaMalloc((void**)&d_q_two[1], sizeof(double)*2*M)); // for next
        gpu_error_check(cudaMalloc((void**)&d_q_half_step, sizeof(double)*M));
        gpu_error_check(cudaMalloc((void**)&d_q_junction,  sizeof(double)*M));

        double q_unity[M];
        for(int i=0; i<M; i++)
            q_unity[i] = 1.0;
        gpu_error_check(cudaMalloc((void**)&d_q_unity, sizeof(double)*M));
        gpu_error_check(cudaMemcpy(d_q_unity, q_unity, sizeof(double)*M, cudaMemcpyHostToDevice));

        // Allocate memory for stress calculation: compute_stress()
        for(int gpu=0; gpu<N_GPUS; gpu++)
        {
            gpu_error_check(cudaSetDevice(gpu));
            gpu_error_check(cudaMalloc((void**)&d_stress_q[gpu][0],     sizeof(double)*2*M)); // prev
            gpu_error_check(cudaMalloc((void**)&d_stress_q[gpu][1],     sizeof(double)*2*M)); // next
        }

        gpu_error_check(cudaSetDevice(0));
        // For concentration computation
        gpu_error_check(cudaMalloc((void**)&d_q_block_v[0], sizeof(double)*M)); // for prev
        gpu_error_check(cudaMalloc((void**)&d_q_block_v[1], sizeof(double)*M)); // for next
        gpu_error_check(cudaMalloc((void**)&d_q_block_u[0], sizeof(double)*M)); // for prev
        gpu_error_check(cudaMalloc((void**)&d_q_block_u[1], sizeof(double)*M)); // for next
        gpu_error_check(cudaMalloc((void**)&d_phi,          sizeof(double)*M));

        propagator_solver->update_bond_function();
        gpu_error_check(cudaSetDevice(0));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
CudaSolverReduceMemoryDiscrete::~CudaSolverReduceMemoryDiscrete()
{
    const int N_GPUS = CudaCommon::get_instance().get_n_gpus();
    
    delete sc;

    delete[] single_polymer_partitions;
    delete[] single_solvent_partitions;

    for(const auto& item: propagator)
    {
        for(int i=0; i<propagator_size[item.first]; i++)
            cudaFreeHost(item.second[i]);
        delete[] item.second;
    }
    for(const auto& item: phi_block)
        cudaFreeHost(item.second);
    for(const auto& item: propagator_junction)
        cudaFreeHost(item.second);
    for(const auto& item: phi_solvent)
        delete[] item;

    #ifndef NDEBUG
    for(const auto& item: propagator_finished)
        delete[] item.second;
    #endif

    for(int gpu=0; gpu<N_GPUS; gpu++)
    {
        cudaFree(d_q_one[gpu][0]);
        cudaFree(d_q_one[gpu][1]);
        cudaFree(d_q_mask[gpu]);
    }
    cudaFree(d_q_two[0]);
    cudaFree(d_q_two[1]);
    cudaFree(d_phi);
    cudaFree(d_q_unity);
    cudaFree(d_q_half_step);
    cudaFree(d_q_junction);

    // For stress calculation: compute_stress()
    for(int gpu=0; gpu<N_GPUS; gpu++)
    {
        cudaFree(d_stress_q[gpu][0]);
        cudaFree(d_stress_q[gpu][1]);
    }

    // For concentration computation
    cudaFree(d_q_block_v[0]);
    cudaFree(d_q_block_v[1]);
    cudaFree(d_q_block_u[0]);
    cudaFree(d_q_block_u[1]);

    // Destroy streams
    for(int gpu=0; gpu<N_GPUS; gpu++)
    {
        cudaStreamDestroy(streams[gpu][0]);
        cudaStreamDestroy(streams[gpu][1]);
    }
}

void CudaSolverReduceMemoryDiscrete::update_bond_function()
{
    try
    {
        propagator_solver->update_bond_function();
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaSolverReduceMemoryDiscrete::compute_statistics(
    std::string device,
    std::map<std::string, const double*> w_input,
    std::map<std::string, const double*> q_init)
{
    try
    {
        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();
        const int N_GPUS = CudaCommon::get_instance().get_n_gpus();

        const int M = cb->get_n_grid();
        const double ds = molecules->get_ds();

        cudaMemcpyKind cudaMemcpyInputToDevice;
        if (device == "gpu")
            cudaMemcpyInputToDevice = cudaMemcpyDeviceToDevice;
        else if(device == "cpu")
            cudaMemcpyInputToDevice = cudaMemcpyHostToDevice;
        else
        {
            throw_with_line_number("Invalid device \"" + device + "\".");
        }

        for(const auto& item: propagators_analyzer->get_essential_propagator_codes())
        {
            if( w_input.find(item.second.monomer_type) == w_input.end())
                throw_with_line_number("monomer_type \"" + item.second.monomer_type + "\" is not in w_input.");
        }

        // Copy mask to d_q_mask
        for(int gpu=0; gpu<N_GPUS; gpu++)
        {
            gpu_error_check(cudaSetDevice(gpu));
            if (cb->get_mask() != nullptr)
            {
                gpu_error_check(cudaMemcpy(d_q_mask[gpu], cb->get_mask(), sizeof(double)*M, cudaMemcpyInputToDevice));
            }
            else
            {
                d_q_mask[gpu] = nullptr;
            }
        }

        // Update dw or d_exp_dw
        propagator_solver->update_dw(device, w_input);

        // For each time span
        auto& branch_schedule = sc->get_schedule();
        for (auto parallel_job = branch_schedule.begin(); parallel_job != branch_schedule.end(); parallel_job++)
        {
            gpu_error_check(cudaSetDevice(0));
            // For each propagator
            for(size_t job=0; job<parallel_job->size(); job++)
            {
                auto& key = std::get<0>((*parallel_job)[job]);
                int n_segment_from = std::get<1>((*parallel_job)[job]);
                int n_segment_to = std::get<2>((*parallel_job)[job]);
                auto& deps = propagators_analyzer->get_essential_propagator_code(key).deps;
                auto monomer_type = propagators_analyzer->get_essential_propagator_code(key).monomer_type;

                // Check key
                #ifndef NDEBUG
                if (propagator.find(key) == propagator.end())
                    throw_with_line_number("Could not find key '" + key + "'. ");
                #endif

                double **_propagator = propagator[key];
                double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

                // Calculate one block end
                if(n_segment_from == 1 && deps.size() == 0) // if it is leaf node
                {
                    // q_init
                    if (key[0] == '{')
                    {
                        std::string g = PropagatorCode::get_q_input_idx_from_key(key);
                        if (q_init.find(g) == q_init.end())
                            throw_with_line_number( "Could not find q_init[\"" + g + "\"].");
                        gpu_error_check(cudaMemcpy(d_q_one[0][0], q_init[g], sizeof(double)*M, cudaMemcpyInputToDevice));
                        multi_real<<<N_BLOCKS, N_THREADS>>>(d_q_one[0][0], d_q_one[0][0], _d_exp_dw, 1.0, M);
                    }
                    else
                    {
                        gpu_error_check(cudaMemcpy(d_q_one[0][0], _d_exp_dw, sizeof(double)*M, cudaMemcpyDeviceToDevice));
                    }

                    // Multiply mask
                    if (d_q_mask[0] != nullptr)
                        multi_real<<<N_BLOCKS, N_THREADS>>>(d_q_one[0][0], d_q_one[0][0], d_q_mask[0], 1.0, M);

                    gpu_error_check(cudaMemcpy(_propagator[0], d_q_one[0][0], sizeof(double)*M, cudaMemcpyDeviceToHost));

                    #ifndef NDEBUG
                    propagator_finished[key][0] = true;
                    #endif
                }
                else if (n_segment_from == 1 && deps.size() > 0) // if it is not leaf node
                {
                    // If it is aggregated
                    if (key[0] == '[')
                    {
                        // Initialize to zero
                        gpu_error_check(cudaMemset(d_q_one[0][0], 0, sizeof(double)*M));

                        int prev, next;
                        prev = 0;
                        next = 1;

                        // Copy memory from host to device
                        std::string sub_dep = std::get<0>(deps[0]);
                        int sub_n_segment   = std::get<1>(deps[0]);
                        int sub_n_repeated;
                        gpu_error_check(cudaMemcpy(d_propagator_sub_dep[prev], propagator[sub_dep][sub_n_segment-1], sizeof(double)*M, cudaMemcpyHostToDevice));

                        for(size_t d=0; d<deps.size(); d++)
                        {
                            sub_dep         = std::get<0>(deps[d]);
                            sub_n_segment   = std::get<1>(deps[d]);
                            sub_n_repeated  = std::get<2>(deps[d]);

                            // Check sub key
                            #ifndef NDEBUG
                            if (propagator.find(sub_dep) == propagator.end())
                                throw_with_line_number("Could not find sub key '" + sub_dep + "'. ");
                            if (!propagator_finished[sub_dep][sub_n_segment-1])
                                throw_with_line_number("Could not compute '" + key +  "', since '"+ sub_dep + std::to_string(sub_n_segment) + "' is not prepared.");
                            #endif

                            // STREAM 1: copy memory from host to device
                            if (d < deps.size()-1)
                            {
                                std::string sub_dep_next = std::get<0>(deps[d+1]);
                                int sub_n_segment_next   = std::get<1>(deps[d+1]);

                                gpu_error_check(cudaMemcpyAsync(d_propagator_sub_dep[next],
                                                propagator[sub_dep_next][sub_n_segment_next-1], sizeof(double)*M,
                                                cudaMemcpyHostToDevice, streams[0][1]));
                            }

                            // STREAM 0: compute linear combination
                            lin_comb<<<N_BLOCKS, N_THREADS, 0, streams[0][0]>>>(
                                    d_q_one[0][0], 1.0, d_q_one[0][0],
                                    sub_n_repeated, d_propagator_sub_dep[prev], M);

                            std::swap(prev, next);
                            cudaDeviceSynchronize();
                        }
                        propagator_solver->advance_one_propagator_discrete(0,
                            d_q_one[0][0], d_q_one[0][0],
                            monomer_type, d_q_mask[0]);

                        gpu_error_check(cudaMemcpy(_propagator[0], d_q_one[0][0], sizeof(double)*M, cudaMemcpyDeviceToHost));

                        #ifndef NDEBUG
                        propagator_finished[key][0] = true;
                        #endif
                    }
                    else
                    { 
                        // Initialize to one
                        gpu_error_check(cudaMemcpy(d_q_junction, d_q_unity, sizeof(double)*M, cudaMemcpyDeviceToDevice));

                        int prev, next;
                        prev = 0;
                        next = 1;

                        // Copy memory from host to device
                        std::string sub_dep = std::get<0>(deps[0]);
                        int sub_n_segment   = std::get<1>(deps[0]);
                        gpu_error_check(cudaMemcpy(d_propagator_sub_dep[prev], propagator[sub_dep][sub_n_segment-1], sizeof(double)*M, cudaMemcpyHostToDevice));

                        for(size_t d=0; d<deps.size(); d++)
                        {
                            sub_dep       = std::get<0>(deps[d]);
                            sub_n_segment = std::get<1>(deps[d]);

                            // Check sub key
                            #ifndef NDEBUG
                            if (propagator.find(sub_dep) == propagator.end())
                                throw_with_line_number("Could not find sub key '" + sub_dep + "'. ");
                            if (!propagator_finished[sub_dep][sub_n_segment-1])
                                throw_with_line_number("Could not compute '" + key +  "', since '"+ sub_dep + std::to_string(sub_n_segment) + "' is not prepared.");
                            #endif

                            // STREAM 1: copy memory from host to device
                            if (d < deps.size()-1)
                            {
                                std::string sub_dep_next = std::get<0>(deps[d+1]);
                                int sub_n_segment_next   = std::get<1>(deps[d+1]);

                                gpu_error_check(cudaMemcpyAsync(d_propagator_sub_dep[next],
                                                propagator[sub_dep_next][sub_n_segment_next-1], sizeof(double)*M,
                                                cudaMemcpyHostToDevice, streams[0][1]));
                            }

                            // STREAM 0: advance propagator half step
                            propagator_solver->advance_propagator_discrete_half_bond_step(0,
                                d_propagator_sub_dep[prev],
                                d_q_half_step, propagators_analyzer->get_essential_propagator_code(sub_dep).monomer_type);
                            multi_real<<<N_BLOCKS, N_THREADS, 0, streams[0][0]>>>(d_q_junction, d_q_junction, d_q_half_step, 1.0, M);

                            std::swap(prev, next);
                            cudaDeviceSynchronize();
                        }
                        gpu_error_check(cudaMemcpy(propagator_junction[key], d_q_junction, sizeof(double)*M,cudaMemcpyDeviceToHost));

                        // Add half bond
                        propagator_solver->advance_propagator_discrete_half_bond_step(0, d_q_junction, d_q_one[0][0], monomer_type);

                        // Add full segment
                        multi_real<<<N_BLOCKS, N_THREADS>>>(d_q_one[0][0], d_q_one[0][0], _d_exp_dw, 1.0, M);

                        // Multiply mask
                        if (d_q_mask[0] != nullptr)
                            multi_real<<<N_BLOCKS, N_THREADS>>>(d_q_one[0][0], d_q_one[0][0], d_q_mask[0], 1.0, M);

                        gpu_error_check(cudaMemcpy(_propagator[0], d_q_one[0][0], sizeof(double)*M, cudaMemcpyDeviceToHost));

                        #ifndef NDEBUG
                        propagator_finished[key][0] = true;
                        #endif
                    }
                }
                else
                {
                    int n = n_segment_from-1;

                    #ifndef NDEBUG
                    if (!propagator_finished[key][n-1])
                        throw_with_line_number("unfinished, key: " + key + ", " + std::to_string(n-1));
                    #endif
            
                    //  copy memory from host to device 
                    gpu_error_check(cudaMemcpy(d_q_one[0][0], _propagator[n-1], sizeof(double)*M,
                        cudaMemcpyHostToDevice));

                    propagator_solver->advance_one_propagator_discrete(0, 
                        d_q_one[0][0],
                        d_q_one[0][0],
                        monomer_type, d_q_mask[0]);

                    cudaDeviceSynchronize();

                    // Copy memory from device to host
                    gpu_error_check(cudaMemcpy(_propagator[n], d_q_one[0][0], sizeof(double)*M,
                        cudaMemcpyDeviceToHost));

                    #ifndef NDEBUG
                    propagator_finished[key][n] = true;
                    #endif
                }
                cudaDeviceSynchronize();
            }
            // Synchronize all GPUs
            for(int gpu=0; gpu<N_GPUS; gpu++)
            {
                gpu_error_check(cudaSetDevice(gpu));
                gpu_error_check(cudaDeviceSynchronize());
            }

            // Advance propagator successively
            if(N_GPUS > 1 && parallel_job->size() == 2)
            {
                const int N_JOBS = parallel_job->size();
                std::string keys[N_JOBS];
                int n_segment_froms[N_JOBS];
                int n_segment_tos[N_JOBS];
                std::string monomer_types[N_JOBS];
                double **_propagator_keys[N_JOBS];
                
                for(int j=0; j<N_JOBS; j++)
                {
                    keys[j] = std::get<0>((*parallel_job)[j]);
                    n_segment_froms[j] = std::get<1>((*parallel_job)[j]);
                    n_segment_tos[j] = std::get<2>((*parallel_job)[j]);
                    monomer_types[j] = propagators_analyzer->get_essential_propagator_code(keys[j]).monomer_type;
                    _propagator_keys[j] = propagator[keys[j]];
                }

                int prev, next;
                prev = 0;
                next = 1;

                // Copy propagators from host to device
                for(int gpu=0; gpu<N_GPUS; gpu++)
                    gpu_error_check(cudaMemcpy(d_q_one[gpu][prev], _propagator_keys[gpu][n_segment_froms[gpu]-1], sizeof(double)*M,
                        cudaMemcpyHostToDevice));

                for(int n=0; n<n_segment_tos[0]-n_segment_froms[0]; n++)
                {
                    #ifndef NDEBUG
                    for(int gpu=0; gpu<N_GPUS; gpu++)
                    {
                        if (!propagator_finished[keys[gpu]][n-1+n_segment_froms[gpu]])
                            throw_with_line_number("unfinished, key: " + keys[gpu] + ", " + std::to_string(n-1+n_segment_froms[gpu]));
                    }
                    #endif

                    // DEVICE 0,1, STREAM 0: calculate propagators 
                    propagator_solver->advance_two_propagators_discrete_two_gpus(
                        d_q_one[0][prev], d_q_one[1][prev],
                        d_q_one[0][next], d_q_one[1][next],
                        monomer_types[0], monomer_types[1],
                        d_q_mask);

                    // STREAM 1: copy propagators from device to host
                    for(int gpu=0; gpu<N_GPUS; gpu++)
                    {
                        gpu_error_check(cudaSetDevice(gpu));
                        if (n > 0)
                        {
                            gpu_error_check(cudaMemcpyAsync(_propagator_keys[gpu][n-1+n_segment_froms[gpu]], d_q_one[gpu][prev], sizeof(double)*M,
                                cudaMemcpyDeviceToHost, streams[gpu][1]));
                        }
                    }

                    // Synchronize all GPUs
                    for(int gpu=0; gpu<N_GPUS; gpu++)
                    {
                        gpu_error_check(cudaSetDevice(gpu));
                        gpu_error_check(cudaDeviceSynchronize());
                    }
                    std::swap(prev, next);

                    #ifndef NDEBUG
                    for(int gpu=0; gpu<N_GPUS; gpu++)
                        propagator_finished[keys[gpu]][n+n_segment_froms[gpu]] = true;
                    #endif
                }
                for(int gpu=0; gpu<N_GPUS; gpu++)
                {
                    // Copy propagators from device to host
                    gpu_error_check(cudaMemcpy(_propagator_keys[gpu][n_segment_tos[gpu]-1], d_q_one[gpu][prev], sizeof(double)*M,
                        cudaMemcpyDeviceToHost));
                }
            }
            else if (N_GPUS == 1 && parallel_job->size() == 2)
            {
                const int N_JOBS = parallel_job->size();
                std::string keys[N_JOBS];
                int n_segment_froms[N_JOBS];
                int n_segment_tos[N_JOBS];
                std::string monomer_types[N_JOBS];
                double **_propagator_keys[N_JOBS];
                
                for(int j=0; j<N_JOBS; j++)
                {
                    keys[j] = std::get<0>((*parallel_job)[j]);
                    n_segment_froms[j] = std::get<1>((*parallel_job)[j]);
                    n_segment_tos[j] = std::get<2>((*parallel_job)[j]);
                    monomer_types[j] = propagators_analyzer->get_essential_propagator_code(keys[j]).monomer_type;
                    _propagator_keys[j] = propagator[keys[j]];
                }

                int prev, next;
                prev = 0;
                next = 1;

                gpu_error_check(cudaSetDevice(0));
                //  copy memory from host to device
                for(int j=0; j<N_JOBS; j++)
                    gpu_error_check(cudaMemcpy(&d_q_two[prev][j*M], _propagator_keys[j][n_segment_froms[j]-1], sizeof(double)*M,
                        cudaMemcpyDeviceToHost));

                for(int n=0; n<n_segment_tos[0]-n_segment_froms[0]; n++)
                {
                    #ifndef NDEBUG
                    for(int j=0; j<N_JOBS; j++)
                    {
                        if (!propagator_finished[keys[j]][n-1+n_segment_froms[j]])
                            throw_with_line_number("unfinished, key: " + keys[j] + ", " + std::to_string(n-n_segment_froms[j]));
                    }
                    #endif

                    // STREAM 0: calculate propagator
                    propagator_solver->advance_two_propagators_discrete_without_copy(
                        d_q_two[prev], d_q_two[next],
                        monomer_types[0], monomer_types[1],
                        d_q_mask[0]);

                    // STREAM 1: copy memory from device to host
                    if (n > 0)
                    {
                        for(int j=0; j<N_JOBS; j++)
                            gpu_error_check(cudaMemcpyAsync(_propagator_keys[j][n-1+n_segment_froms[j]], &d_q_two[prev][j*M], sizeof(double)*M,
                                cudaMemcpyHostToDevice, streams[0][1]));
                    }

                    std::swap(prev, next);
                    cudaDeviceSynchronize();

                    #ifndef NDEBUG
                    for(int j=0; j<N_JOBS; j++)
                        propagator_finished[keys[j]][n+n_segment_froms[j]] = true;
                    #endif
                }
                // Copy memory from device to host
                for(int j=0; j<N_JOBS; j++)
                    gpu_error_check(cudaMemcpy(_propagator_keys[j][n_segment_tos[j]-1], &d_q_two[prev][j*M], sizeof(double)*M,
                        cudaMemcpyDeviceToHost));
            }
            else if (parallel_job->size() > 0)
            {
                const int N_JOBS = parallel_job->size();
                std::string keys[N_JOBS];
                int n_segment_froms[N_JOBS];
                int n_segment_tos[N_JOBS];
                std::string monomer_types[N_JOBS];
                double **_propagator_keys[N_JOBS];
                
                for(int j=0; j<N_JOBS; j++)
                {
                    keys[j] = std::get<0>((*parallel_job)[j]);
                    n_segment_froms[j] = std::get<1>((*parallel_job)[j]);
                    n_segment_tos[j] = std::get<2>((*parallel_job)[j]);
                    monomer_types[j] = propagators_analyzer->get_essential_propagator_code(keys[j]).monomer_type;
                    _propagator_keys[j] = propagator[keys[j]];
                }

                for(int j=0; j<N_JOBS; j++)
                {
                    int prev, next;
                    prev = 0;
                    next = 1;

                    gpu_error_check(cudaSetDevice(0));
                    //  copy memory from host to device 
                    gpu_error_check(cudaMemcpy(d_q_one[0][prev], _propagator_keys[j][n_segment_froms[j]-1], sizeof(double)*M,
                        cudaMemcpyHostToDevice));

                    for(int n=n_segment_froms[j]; n<n_segment_tos[j]; n++)
                    {
                        #ifndef NDEBUG
                        if (!propagator_finished[keys[j]][n-1])
                            throw_with_line_number("unfinished, key: " + keys[j] + ", " + std::to_string(n));
                        #endif

                        // STREAM 1: copy memory from device to host
                        if (n > n_segment_froms[j])
                        {
                            gpu_error_check(cudaMemcpyAsync(_propagator_keys[j][n-1], d_q_one[0][prev], sizeof(double)*M,
                                cudaMemcpyDeviceToHost, streams[0][1]));
                        }

                        // STREAM 0: calculate propagator
                        propagator_solver->advance_one_propagator_discrete(0,
                            d_q_one[0][prev],
                            d_q_one[0][next],
                            monomer_types[j], d_q_mask[0]);

                        std::swap(prev, next);
                        cudaDeviceSynchronize();

                        #ifndef NDEBUG
                        propagator_finished[keys[j]][n] = true;
                        #endif
                    }
                    // Copy memory from device to host
                    gpu_error_check(cudaMemcpy(_propagator_keys[j][n_segment_tos[j]-1], d_q_one[0][prev], sizeof(double)*M,
                        cudaMemcpyDeviceToHost));
                }
            }
            // Synchronize all GPUs
            for(int gpu=0; gpu<N_GPUS; gpu++)
            {
                gpu_error_check(cudaSetDevice(gpu));
                gpu_error_check(cudaDeviceSynchronize());
            }
        }
        gpu_error_check(cudaSetDevice(0));

        // Compute total partition function of each distinct polymers
        for(const auto& segment_info: single_partition_segment)
        {
            int p                  = std::get<0>(segment_info);
            double *propagator_v   = std::get<1>(segment_info);
            double *propagator_u   = std::get<2>(segment_info);
            std::string monomer_type = std::get<3>(segment_info);
            int n_aggregated         = std::get<4>(segment_info);
            double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

            // Copy propagators from host to device
            gpu_error_check(cudaMemcpy(d_q_block_v[0], propagator_v, sizeof(double)*M, cudaMemcpyHostToDevice));
            gpu_error_check(cudaMemcpy(d_q_block_u[0], propagator_u, sizeof(double)*M, cudaMemcpyHostToDevice));

            single_polymer_partitions[p] = cb->inner_product_inverse_weight_device(
                d_q_block_v[0],  // q
                d_q_block_u[0],  // q^dagger
                _d_exp_dw)/n_aggregated/cb->get_volume();
        }

        // Calculate segment concentrations
        for(const auto& block: phi_block)
        {
            const auto& key = block.first;
            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            int n_repeated;
            int n_segment_allocated = propagators_analyzer->get_essential_block(key).n_segment_allocated;
            int n_segment_offset    = propagators_analyzer->get_essential_block(key).n_segment_offset;
            int n_segment_original  = propagators_analyzer->get_essential_block(key).n_segment_original;
            std::string monomer_type = propagators_analyzer->get_essential_block(key).monomer_type;
            double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

            // Contains no '['
            if (dep_u.find('[') == std::string::npos)
                n_repeated = propagators_analyzer->get_essential_block(key).v_u.size();
            else
                n_repeated = 1;

            // Check keys
            #ifndef NDEBUG
            if (propagator.find(dep_v) == propagator.end())
                throw_with_line_number("Could not find dep_v key'" + dep_v + "'. ");
            if (propagator.find(dep_u) == propagator.end())
                throw_with_line_number("Could not find dep_u key'" + dep_u + "'. ");
            #endif

            // Normalize concentration
            Polymer& pc = molecules->get_polymer(p);
            double norm = molecules->get_ds()*pc.get_volume_fraction()/pc.get_alpha()/single_polymer_partitions[p]*n_repeated;

            // Calculate phi of one block (possibly multiple blocks when using aggregation)
            calculate_phi_one_block(
                block.second,               // phi
                propagator[dep_v],          // dependency v
                propagator[dep_u],          // dependency u
                _d_exp_dw,                  // exp_dw
                n_segment_allocated,
                n_segment_offset,
                n_segment_original,
                norm);
        }
        // Calculate partition functions and concentrations of solvents
        for(size_t s=0; s<molecules->get_n_solvent_types(); s++)
        {
            double volume_fraction = std::get<0>(molecules->get_solvent(s));
            std::string monomer_type = std::get<1>(molecules->get_solvent(s));
            double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

            single_solvent_partitions[s] = cb->integral_device(_d_exp_dw)/cb->get_volume();
            linear_scaling_real<<<N_BLOCKS, N_THREADS>>>(d_phi, _d_exp_dw, volume_fraction/single_solvent_partitions[s], 0.0, M);
            gpu_error_check(cudaMemcpy(phi_solvent[s], d_phi, sizeof(double)*M, cudaMemcpyDeviceToHost));
        }
        gpu_error_check(cudaSetDevice(0));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaSolverReduceMemoryDiscrete::calculate_phi_one_block(
    double *phi, double **q_1, double **q_2, double *d_exp_dw, const int N, const int N_OFFSET, const int N_ORIGINAL, const double NORM)
{
    try
    {
        gpu_error_check(cudaSetDevice(0));
        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();
        const int M = cb->get_n_grid();

        int prev, next;
        prev = 0;
        next = 1;

        // Copy propagators from host to device
        gpu_error_check(cudaMemcpy(d_q_block_v[prev], q_1[N_ORIGINAL-N_OFFSET-1], sizeof(double)*M, cudaMemcpyHostToDevice));
        gpu_error_check(cudaMemcpy(d_q_block_u[prev], q_2[0],                     sizeof(double)*M, cudaMemcpyHostToDevice));

        // Initialize to zero
        gpu_error_check(cudaMemset(d_phi, 0, sizeof(double)*M));
 
        for(int n=0; n<N; n++)
        {
            // STREAM 1: copy propagators from host to device
            if (n+1 < N)
            {
                gpu_error_check(cudaMemcpyAsync(d_q_block_v[next], q_1[N_ORIGINAL-N_OFFSET-(n+1)-1],
                    sizeof(double)*M, cudaMemcpyHostToDevice, streams[0][1]));
                gpu_error_check(cudaMemcpyAsync(d_q_block_u[next], q_2[(n+1)],
                    sizeof(double)*M, cudaMemcpyHostToDevice, streams[0][1]));
            }

            // STREAM 0: multiply two propagators
            add_multi_real<<<N_BLOCKS, N_THREADS, 0, streams[0][0]>>>(d_phi, d_q_block_v[prev], d_q_block_u[prev], NORM, M);
            std::swap(prev, next);
            cudaDeviceSynchronize();
        }
        // divide exp_dw
        divide_real<<<N_BLOCKS, N_THREADS>>>(d_phi, d_phi, d_exp_dw, 1.0, M);

        // Copy propagators from device to host
        gpu_error_check(cudaMemcpy(phi, d_phi, sizeof(double)*M, cudaMemcpyDeviceToHost));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
double CudaSolverReduceMemoryDiscrete::get_total_partition(int polymer)
{
    try
    {
        return single_polymer_partitions[polymer];
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaSolverReduceMemoryDiscrete::get_total_concentration(std::string monomer_type, double *phi)
{
    try
    {
        const int M = cb->get_n_grid();
        // Initialize array
        for(int i=0; i<M; i++)
            phi[i] = 0.0;

        // For each block
        for(const auto& block: phi_block)
        {
            std::string dep_v = std::get<1>(block.first);
            int n_segment_allocated = propagators_analyzer->get_essential_block(block.first).n_segment_allocated;
            if (PropagatorCode::get_monomer_type_from_key(dep_v) == monomer_type && n_segment_allocated != 0)
            {
                for(int i=0; i<M; i++)
                    phi[i] += block.second[i]; 
            }
        }
        // For each solvent
        for(int s=0;s<molecules->get_n_solvent_types();s++)
        {
            if (std::get<1>(molecules->get_solvent(s)) == monomer_type)
            {
                double *phi_solvent_ = phi_solvent[s];
                for(int i=0; i<M; i++)
                    phi[i] += phi_solvent_[i];
            }
        }
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaSolverReduceMemoryDiscrete::get_total_concentration(int p, std::string monomer_type, double *phi)
{
    try
    {
        const int M = cb->get_n_grid();
        const int P = molecules->get_n_polymer_types();

        if (p < 0 || p > P-1)
            throw_with_line_number("Index (" + std::to_string(p) + ") must be in range [0, " + std::to_string(P-1) + "]");

        // Initialize array
        for(int i=0; i<M; i++)
            phi[i] = 0.0;

        // For each block
        for(const auto& block: phi_block)
        {
            int polymer_idx = std::get<0>(block.first);
            std::string dep_v = std::get<1>(block.first);
            int n_segment_allocated = propagators_analyzer->get_essential_block(block.first).n_segment_allocated;
            if (polymer_idx == p && PropagatorCode::get_monomer_type_from_key(dep_v) == monomer_type && n_segment_allocated != 0)
            {
                for(int i=0; i<M; i++)
                    phi[i] += block.second[i]; 
            }
        }
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaSolverReduceMemoryDiscrete::get_block_concentration(int p, double *phi)
{
    try
    {
        const int M = cb->get_n_grid();
        const int P = molecules->get_n_polymer_types();

        if (p < 0 || p > P-1)
            throw_with_line_number("Index (" + std::to_string(p) + ") must be in range [0, " + std::to_string(P-1) + "]");

        if (propagators_analyzer->is_using_propagator_aggregation())
            throw_with_line_number("Disable 'aggregation' option to obtain concentration of each block.");

        Polymer& pc = molecules->get_polymer(p);
        std::vector<Block>& blocks = pc.get_blocks();

        for(size_t b=0; b<blocks.size(); b++)
        {
            std::string dep_v = pc.get_propagator_key(blocks[b].v, blocks[b].u);
            std::string dep_u = pc.get_propagator_key(blocks[b].u, blocks[b].v);
            if (dep_v < dep_u)
                dep_v.swap(dep_u);

            double* _essential_phi_block = phi_block[std::make_tuple(p, dep_v, dep_u)];
            for(int i=0; i<M; i++)
                phi[i+b*M] = _essential_phi_block[i]; 
        }
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
double CudaSolverReduceMemoryDiscrete::get_solvent_partition(int s)
{
    try
    {
        return single_solvent_partitions[s];
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaSolverReduceMemoryDiscrete::get_solvent_concentration(int s, double *phi)
{
    try
    {
        const int M = cb->get_n_grid();
        const int S = molecules->get_n_solvent_types();

        if (s < 0 || s > S-1)
            throw_with_line_number("Index (" + std::to_string(s) + ") must be in range [0, " + std::to_string(S-1) + "]");

        double *phi_solvent_ = phi_solvent[s];
        for(int i=0; i<M; i++)
            phi[i] = phi_solvent_[i];
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
std::vector<double> CudaSolverReduceMemoryDiscrete::compute_stress()
{
    // This method should be invoked after invoking compute_statistics().

    // To calculate stress, we multiply weighted fourier basis to q(k)*q^dagger(-k).
    // We only need the real part of stress calculation.
    try
    {
        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();
        const int N_GPUS = CudaCommon::get_instance().get_n_gpus();

        const int DIM  = cb->get_dim();
        const int M    = cb->get_n_grid();

        std::vector<double> stress(DIM);
        std::map<std::tuple<int, std::string, std::string>, std::array<double,3>> block_dq_dl[MAX_GPUS];

        // Compute stress for each block
        for(const auto& block: phi_block)
        {
            const auto& key = block.first;
            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            const int N           = propagators_analyzer->get_essential_block(key).n_segment_allocated;
            const int N_OFFSET    = propagators_analyzer->get_essential_block(key).n_segment_offset;
            const int N_ORIGINAL  = propagators_analyzer->get_essential_block(key).n_segment_original;
            std::string monomer_type = propagators_analyzer->get_essential_block(key).monomer_type;

            // Contains no '['
            int n_repeated;
            if (dep_u.find('[') == std::string::npos)
                n_repeated = propagators_analyzer->get_essential_block(key).v_u.size();
            else
                n_repeated = 1;

            double **q_1 = propagator[dep_v];         // Propagator q
            double **q_2 = propagator[dep_u];         // Propagator q^dagger

            std::array<double,3> _block_dq_dl[MAX_GPUS];
            for(int gpu=0; gpu<N_GPUS; gpu++)
                for(int d=0; d<3; d++)
                    _block_dq_dl[gpu][d] = 0.0;

            // Check block_stress_info
            const auto& _block_stress_info_key = block_stress_info[key];
            if(_block_stress_info_key.size() != (unsigned int) (N+1))
            {
                throw_with_line_number("Mismatch of block_stress_info("
                    + std::to_string(p) + "," + dep_v + "," + dep_u + ") "
                    + std::to_string(_block_stress_info_key.size()) + ".size() with N+1 (" + std::to_string(N+1) + ")");
            }

            // Variables for block_stress_info
            double *propagator_v;
            double *propagator_u;

            int prev, next;
            prev = 0;
            next = 1;

            // Copy memory from device to device
            for(int gpu=0; gpu<N_GPUS; gpu++)
            {
                // Index
                int idx = gpu;
                gpu_error_check(cudaSetDevice(gpu));
                if (idx <= N)
                {
                    propagator_v = std::get<0>(_block_stress_info_key[idx]);
                    propagator_u = std::get<1>(_block_stress_info_key[idx]);

                    if (propagator_v != nullptr)
                    {
                        gpu_error_check(cudaMemcpy(&d_stress_q[gpu][prev][0], propagator_v, sizeof(double)*M, cudaMemcpyHostToDevice));
                        gpu_error_check(cudaMemcpy(&d_stress_q[gpu][prev][M], propagator_u, sizeof(double)*M, cudaMemcpyHostToDevice));
                    }
                }
            }

            // Compute stress
            for(int n=0; n<=N; n+=N_GPUS)
            {
                // STREAM 1: copy memory from device to device
                for(int gpu=0; gpu<N_GPUS; gpu++)
                {
                    const int idx = n + gpu;
                    const int idx_next = idx + N_GPUS;
                    gpu_error_check(cudaSetDevice(gpu));
                    if (idx_next <= N)
                    {
                        propagator_v = std::get<0>(_block_stress_info_key[idx_next]);
                        propagator_u = std::get<1>(_block_stress_info_key[idx_next]);

                        if (propagator_v != nullptr)
                        {
                            gpu_error_check(cudaMemcpyAsync(&d_stress_q[gpu][next][0], propagator_v, sizeof(double)*M, cudaMemcpyHostToDevice, streams[gpu][1]));
                            gpu_error_check(cudaMemcpyAsync(&d_stress_q[gpu][next][M], propagator_u, sizeof(double)*M, cudaMemcpyHostToDevice, streams[gpu][1]));
                        }
                    }
                }

                // STREAM 0: execute kernels
                for(int gpu=0; gpu<N_GPUS; gpu++)
                {
                    const int idx = n + gpu;
                    gpu_error_check(cudaSetDevice(gpu));
                    if (idx <= N)
                    {
                        propagator_v = std::get<0>(_block_stress_info_key[idx]);

                        if (propagator_v != nullptr)
                        {
                            propagator_solver->compute_single_segment_stress_fourier(gpu, d_stress_q[gpu][prev]);
                        }
                    }
                }

                // STREAM 0: reduction sum
                for(int gpu=0; gpu<N_GPUS; gpu++)
                {
                    const int idx = n + gpu;
                    gpu_error_check(cudaSetDevice(gpu));
                    if (idx <= N)
                    {
                        propagator_v = std::get<0>(_block_stress_info_key[idx]);
                        bool is_half_bond_length = std::get<2>(_block_stress_info_key[idx]);
                        
                        if (propagator_v != nullptr)
                        {
                            // Compute
                            std::vector<double> segment_stress = propagator_solver->compute_single_segment_stress_discrete(gpu, monomer_type, is_half_bond_length);
                            for(int d=0; d<DIM; d++)
                                _block_dq_dl[gpu][d] += segment_stress[d]*n_repeated;
                            
                            // std::cout << "n: " << n << ", " << is_half_bond_length << ", " << segment_stress[0] << std::endl;
                        }
                    }
                }
                // Synchronize all GPUs
                for(int gpu=0; gpu<N_GPUS; gpu++)
                {
                    gpu_error_check(cudaSetDevice(gpu));
                    gpu_error_check(cudaDeviceSynchronize());
                }
                std::swap(prev, next);
            }
            // Copy stress data
            for(int gpu=0; gpu<N_GPUS; gpu++)
                block_dq_dl[gpu][key] = _block_dq_dl[gpu];
        }
        gpu_error_check(cudaSetDevice(0));

        // Compute total stress
        for(int d=0; d<DIM; d++)
            stress[d] = 0.0;
        for(const auto& block: phi_block)
        {
            const auto& key = block.first;
            int p             = std::get<0>(key);
            std::string dep_v = std::get<1>(key);
            std::string dep_u = std::get<2>(key);
            Polymer& pc  = molecules->get_polymer(p);

            for(int gpu=0; gpu<N_GPUS; gpu++)
                for(int d=0; d<DIM; d++)
                    stress[d] += block_dq_dl[gpu][key][d]*pc.get_volume_fraction()/pc.get_alpha()/single_polymer_partitions[p];
        }
        for(int d=0; d<DIM; d++)
            stress[d] /= -3.0*cb->get_lx(d)*M*M/molecules->get_ds();
            
        return stress;
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaSolverReduceMemoryDiscrete::get_chain_propagator(double *q_out, int polymer, int v, int u, int n)
{ 
    // This method should be invoked after invoking compute_statistics()

    // Get chain propagator for a selected polymer, block and direction.
    // This is made for debugging and testing.
    try
    {
        const int M = cb->get_n_grid();
        Polymer& pc = molecules->get_polymer(polymer);
        std::string dep = pc.get_propagator_key(v,u);

        if (propagators_analyzer->get_essential_propagator_codes().find(dep) == propagators_analyzer->get_essential_propagator_codes().end())
            throw_with_line_number("Could not find the propagator code '" + dep + "'. Disable 'aggregation' option to obtain propagators_analyzer.");
            
        const int N = propagators_analyzer->get_essential_propagator_codes()[dep].max_n_segment;
        if (n < 1 || n > N)
            throw_with_line_number("n (" + std::to_string(n) + ") must be in range [1, " + std::to_string(N) + "]");

        double* _propagator = propagator[dep][n-1];
        for(int i=0; i<M; i++)
            q_out[i] = _propagator[i];
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
bool CudaSolverReduceMemoryDiscrete::check_total_partition()
{
    const int M = cb->get_n_grid();
    int n_polymer_types = molecules->get_n_polymer_types();
    std::vector<std::vector<double>> total_partitions;
    for(int p=0;p<n_polymer_types;p++)
    {
        std::vector<double> total_partitions_p;
        total_partitions.push_back(total_partitions_p);
    }

    gpu_error_check(cudaSetDevice(0));
    for(const auto& block: phi_block)
    {
        const auto& key = block.first;
        int p                = std::get<0>(key);
        std::string dep_v    = std::get<1>(key);
        std::string dep_u    = std::get<2>(key);

        int n_aggregated;
        int n_segment_allocated = propagators_analyzer->get_essential_block(key).n_segment_allocated;
        int n_segment_offset    = propagators_analyzer->get_essential_block(key).n_segment_offset;
        int n_segment_original  = propagators_analyzer->get_essential_block(key).n_segment_original;
        std::string monomer_type = propagators_analyzer->get_essential_block(key).monomer_type;
        double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

        // std::cout<< p << ", " << dep_v << ", " << dep_u << ": " << n_segment_original << ", " << n_segment_offset << ", " << n_segment_allocated << std::endl;

        // Contains no '['
        if (dep_u.find('[') == std::string::npos)
            n_aggregated = 1;
        else
            n_aggregated = propagators_analyzer->get_essential_block(key).v_u.size();

        for(int n=0;n<n_segment_allocated;n++)
        {
            // Copy propagators from host to device
            gpu_error_check(cudaMemcpy(d_q_block_v[0], propagator[dep_v][n_segment_original-n_segment_offset-n-1], sizeof(double)*M, cudaMemcpyHostToDevice));
            gpu_error_check(cudaMemcpy(d_q_block_u[0], propagator[dep_u][n], sizeof(double)*M, cudaMemcpyHostToDevice));

            double total_partition = cb->inner_product_inverse_weight_device(
                d_q_block_v[0],  // q
                d_q_block_u[0],  // q^dagger
                _d_exp_dw)/n_aggregated/cb->get_volume();

            // std::cout<< p << ", " << n << ": " << total_partition << std::endl;
            total_partitions[p].push_back(total_partition);
        }
    }

    // Find minimum and maximum of total_partitions
    std::cout<< "Polymer id: maximum,  minimum, and difference of total partitions" << std::endl;
    for(size_t p=0;p<total_partitions.size();p++)
    {
        double max_partition = -1e20;
        double min_partition =  1e20;
        for(size_t n=0;n<total_partitions[p].size();n++)
        {
            if (total_partitions[p][n] > max_partition)
                max_partition = total_partitions[p][n];
            if (total_partitions[p][n] < min_partition)
                min_partition = total_partitions[p][n];
        }
        double diff_partition = abs(max_partition - min_partition);

        std::cout<< "\t" << p << ": " << max_partition << ", " << min_partition << ", " << diff_partition << std::endl;
        if (diff_partition > 1e-7)
            return false;
    }
    return true;
}