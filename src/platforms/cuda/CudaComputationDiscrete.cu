#include <complex>
#include <thrust/reduce.h>
#include <omp.h>
#include <iostream>
#include "CudaComputationDiscrete.h"
#include "CudaComputationBox.h"

CudaComputationDiscrete::CudaComputationDiscrete(
    ComputationBox *cb,
    Molecules *molecules,
    PropagatorAnalyzer *propagator_analyzer)
    : PropagatorComputation(cb, molecules, propagator_analyzer)
{
    try
    {
        const int M = cb->get_n_grid();
        const int N_GPUS = CudaCommon::get_instance().get_n_gpus();

        // Copy streams
        for(int i=0; i<N_STREAMS; i++)
        {
            if (N_GPUS == 1)
            {
                gpu_error_check(cudaSetDevice(0));
            }
            else
            {
                gpu_error_check(cudaSetDevice(i));
            }
            gpu_error_check(cudaStreamCreate(&streams[i][0])); // for kernel execution
            gpu_error_check(cudaStreamCreate(&streams[i][1])); // for memcpy
        }

        this->propagator_solver = new CudaSolverPseudo(cb, molecules, streams, false);

        // Allocate memory for propagators
        gpu_error_check(cudaSetDevice(0));
        if( propagator_analyzer->get_computation_propagator_codes().size() == 0)
            throw_with_line_number("There is no propagator code. Add polymers first.");
        for(const auto& item: propagator_analyzer->get_computation_propagator_codes())
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
            d_propagator[key] = new double*[max_n_segment];
            for(int i=0; i<propagator_size[key]; i++)
                gpu_error_check(cudaMalloc((void**)&d_propagator[key][i], sizeof(double)*M));

            #ifndef NDEBUG
            propagator_finished[key] = new bool[max_n_segment];
            for(int i=0; i<max_n_segment;i++)
                propagator_finished[key][i] = false;
            #endif
        }

        // Allocate memory for propagator_junction, which contain propagator at junction of discrete chain
        for(const auto& item: propagator_analyzer->get_computation_propagator_codes())
        {
            std::string key = item.first;
            d_propagator_junction[key] = nullptr;
            gpu_error_check(cudaMalloc((void**)&d_propagator_junction[key], sizeof(double)*M));
        }

        // Allocate memory for concentrations
        if( propagator_analyzer->get_computation_blocks().size() == 0)
            throw_with_line_number("There is no block. Add polymers first.");
        for(const auto& item: propagator_analyzer->get_computation_blocks())
        {
            d_phi_block[item.first] = nullptr;
            gpu_error_check(cudaMalloc((void**)&d_phi_block[item.first], sizeof(double)*M));
        }

        // Total partition functions for each polymer
        single_polymer_partitions = new double[molecules->get_n_polymer_types()];

        // Remember one segment for each polymer chain to compute total partition function
        int current_p = 0;
        for(const auto& d_block: d_phi_block)
        {
            const auto& key = d_block.first;
            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            // Skip if already found one segment
            if (p != current_p)
                continue;

            int n_aggregated;
            int n_segment_original  = propagator_analyzer->get_computation_block(key).n_segment_original;
            std::string monomer_type = propagator_analyzer->get_computation_block(key).monomer_type;

            // Contains no '['
            if (dep_u.find('[') == std::string::npos)
                n_aggregated = 1;
            else
                n_aggregated = propagator_analyzer->get_computation_block(key).v_u.size();

            single_partition_segment.push_back(std::make_tuple(
                p,
                d_propagator[dep_v][n_segment_original-1],  // q
                d_propagator[dep_u][0],                                      // q_dagger
                monomer_type,       
                n_aggregated                   // how many propagators are aggregated
                ));
            current_p++;
        }

       // Find propagators and bond length for each segment to prepare stress computation
        for(const auto& block: d_phi_block)
        {
            const auto& key = block.first;
            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            const int N           = propagator_analyzer->get_computation_block(key).n_segment_allocated;
            const int N_ORIGINAL  = propagator_analyzer->get_computation_block(key).n_segment_original;

            double **d_q_1 = d_propagator[dep_v];    // dependency v
            double **d_q_2 = d_propagator[dep_u];    // dependency u

            auto& _block_stress_compuation_key = block_stress_computation_plan[key];

            // Find propagators and bond length
            for(int n=0; n<=N; n++)
            {
                double *d_propagator_v = nullptr;
                double *d_propagator_u = nullptr;
                bool is_half_bond_length = false;

                // At v
                if (n == N_ORIGINAL)
                {
                    if (propagator_analyzer->get_computation_propagator_code(dep_v).deps.size() == 0) // if v is leaf node, skip
                    {
                        _block_stress_compuation_key.push_back(std::make_tuple(d_propagator_v, d_propagator_u, is_half_bond_length));
                        continue;
                    }
                    
                    d_propagator_v = d_propagator_junction[dep_v];
                    d_propagator_u = d_q_2[N-1];
                    is_half_bond_length = true;
                }
                // At u
                else if (n == 0 && dep_u.find('[') == std::string::npos){
                    if (propagator_analyzer->get_computation_propagator_code(dep_u).deps.size() == 0) // if u is leaf node, skip
                    {
                        _block_stress_compuation_key.push_back(std::make_tuple(d_propagator_v, d_propagator_u, is_half_bond_length));
                        continue;
                    }

                    d_propagator_v = d_q_1[N_ORIGINAL-1];
                    d_propagator_u = d_propagator_junction[dep_u];
                    is_half_bond_length = true;
                }
                // At aggregation junction
                else if (n == 0)
                {
                    _block_stress_compuation_key.push_back(std::make_tuple(d_propagator_v, d_propagator_u, is_half_bond_length));
                    continue;
                }
                // Within the blocks
                else
                {
                    d_propagator_v = d_q_1[N_ORIGINAL-n-1];
                    d_propagator_u = d_q_2[n-1];
                    is_half_bond_length = false;
                }
                _block_stress_compuation_key.push_back(std::make_tuple(d_propagator_v, d_propagator_u, is_half_bond_length));
            }
        }

        // Total partition functions for each solvent
        single_solvent_partitions = new double[molecules->get_n_solvent_types()];

        // Concentrations for each solvent
        for(int s=0;s<molecules->get_n_solvent_types();s++)
        {
            double *d_phi_;
            gpu_error_check(cudaMalloc((void**)&d_phi_, sizeof(double)*M));
            d_phi_solvent.push_back(d_phi_);
        }

        // Create scheduler for computation of propagator
        sc = new Scheduler(propagator_analyzer->get_computation_propagator_codes(), N_STREAMS); 

        // Allocate memory for pseudo-spectral: advance_propagator()
        double q_unity[M];
        for(int i=0; i<M; i++)
            q_unity[i] = 1.0;
        for(int gpu=0; gpu<N_GPUS; gpu++)
        {
            gpu_error_check(cudaSetDevice(gpu));
            gpu_error_check(cudaMalloc((void**)&d_q_mask[gpu], sizeof(double)*M));
            gpu_error_check(cudaMalloc((void**)&d_q_unity[gpu], sizeof(double)*M));
            gpu_error_check(cudaMemcpy(d_q_unity[gpu], q_unity, sizeof(double)*M, cudaMemcpyHostToDevice));
        }

        for(int gpu=1; gpu<N_GPUS; gpu++)
        {
            gpu_error_check(cudaSetDevice(gpu));
            gpu_error_check(cudaMalloc((void**)&d_propagator_device[gpu][0], sizeof(double)*M));  // prev
            gpu_error_check(cudaMalloc((void**)&d_propagator_device[gpu][1], sizeof(double)*M));  // next
        }

        for(int i=0; i<N_STREAMS; i++)
        {
            if (N_GPUS == 1)
            {
                gpu_error_check(cudaSetDevice(0));
            }
            else
            {
                gpu_error_check(cudaSetDevice(i));
            }
            gpu_error_check(cudaMalloc((void**)&d_q_half_step[i], sizeof(double)*M));
            gpu_error_check(cudaMalloc((void**)&d_q_junction[i],  sizeof(double)*M));
        }

        gpu_error_check(cudaSetDevice(0));
        gpu_error_check(cudaMalloc((void**)&d_phi, sizeof(double)*M));

        // Allocate memory for stress calculation: compute_stress()
        for(int i=0; i<N_STREAMS; i++)
        {
            if (N_GPUS == 1)
            {
                gpu_error_check(cudaSetDevice(0));
            }
            else
            {
                gpu_error_check(cudaSetDevice(i));
            }
            gpu_error_check(cudaMalloc((void**)&d_q_pair[i][0],     sizeof(double)*2*M)); // prev
            gpu_error_check(cudaMalloc((void**)&d_q_pair[i][1],     sizeof(double)*2*M)); // next
        }
        
        propagator_solver->update_laplacian_operator();
        gpu_error_check(cudaSetDevice(0));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
CudaComputationDiscrete::~CudaComputationDiscrete()
{
    const int N_GPUS = CudaCommon::get_instance().get_n_gpus();
    
    delete propagator_solver;
    delete sc;

    delete[] single_polymer_partitions;
    delete[] single_solvent_partitions;

    for(const auto& item: d_propagator)
    {
        for(int i=0; i<propagator_size[item.first]; i++)
            cudaFree(item.second[i]);
        delete[] item.second;
    }
    for(const auto& item: d_phi_block)
        cudaFree(item.second);
    for(const auto& item: d_propagator_junction)
        cudaFree(item.second);
    for(const auto& item: d_phi_solvent)
        cudaFree(item);

    #ifndef NDEBUG
    for(const auto& item: propagator_finished)
        delete[] item.second;
    #endif

    for(int i=0; i<N_STREAMS; i++)
    {
        cudaFree(d_q_half_step[i]);
        cudaFree(d_q_junction[i]);
    }
    cudaFree(d_phi);

    // For pseudo-spectral: advance_propagator()
    for(int gpu=0; gpu<N_GPUS; gpu++)
    {
        cudaFree(d_q_mask[gpu]);
        cudaFree(d_q_unity[gpu]);
    }
    
    for(int gpu=1; gpu<N_GPUS; gpu++)
    {
        cudaFree(d_propagator_device[gpu][0]);
        cudaFree(d_propagator_device[gpu][1]);
    }

    // For stress calculation: compute_stress()
    for(int i=0; i<N_STREAMS; i++)
    {
        cudaFree(d_q_pair[i][0]);
        cudaFree(d_q_pair[i][1]);
    }

    // Destroy streams
    for(int i=0; i<N_STREAMS; i++)
    {
        cudaStreamDestroy(streams[i][0]);
        cudaStreamDestroy(streams[i][1]);
    }
}

void CudaComputationDiscrete::update_laplacian_operator()
{
    try
    {
        propagator_solver->update_laplacian_operator();
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaComputationDiscrete::compute_statistics(
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

        for(const auto& item: propagator_analyzer->get_computation_propagator_codes())
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
            // For each propagator
            #pragma omp parallel for num_threads(N_STREAMS) 
            for(size_t job=0; job<parallel_job->size(); job++)
            {
                int gpu;
                const int STREAM = omp_get_thread_num();
                if (N_GPUS == 1)
                {
                    gpu = 0;
                }
                else
                {
                    gpu = omp_get_thread_num();
                }
                gpu_error_check(cudaSetDevice(gpu));

                auto& key = std::get<0>((*parallel_job)[job]);
                int n_segment_from = std::get<1>((*parallel_job)[job]);
                int n_segment_to = std::get<2>((*parallel_job)[job]);
                auto& deps = propagator_analyzer->get_computation_propagator_code(key).deps;
                auto monomer_type = propagator_analyzer->get_computation_propagator_code(key).monomer_type;
                bool is_initialized = true;

                // Check key
                #ifndef NDEBUG
                if (d_propagator.find(key) == d_propagator.end())
                    throw_with_line_number("Could not find key '" + key + "'. ");
                #endif
                double **_d_propagator = d_propagator[key];
                double *_d_exp_dw = propagator_solver->d_exp_dw[gpu][monomer_type];

                // Calculate one block end
                if(n_segment_from == 1 && deps.size() == 0) // if it is leaf node
                {
                     // q_init
                    if (key[0] == '{')
                    {
                        std::string g = PropagatorCode::get_q_input_idx_from_key(key);
                        if (q_init.find(g) == q_init.end())
                            throw_with_line_number("Could not find q_init[\"" + g + "\"].");
                        gpu_error_check(cudaMemcpy(_d_propagator[0], q_init[g], sizeof(double)*M, cudaMemcpyInputToDevice));
                        multi_real<<<N_BLOCKS, N_THREADS>>>(_d_propagator[0], _d_propagator[0], _d_exp_dw, 1.0, M);
                    }
                    else
                    {
                        gpu_error_check(cudaMemcpy(_d_propagator[0], _d_exp_dw, sizeof(double)*M, cudaMemcpyDeviceToDevice));
                    }
                    
                    #ifndef NDEBUG
                    propagator_finished[key][0] = true;
                    #endif
                }
                // If it is not leaf node
                else if (n_segment_from == 1 && deps.size() > 0)
                {
                    // If it is aggregated
                    if (key[0] == '[')
                    {
                        // Initialize to zero
                        gpu_error_check(cudaMemset(_d_propagator[0], 0, sizeof(double)*M));

                        for(size_t d=0; d<deps.size(); d++)
                        {
                            std::string sub_dep = std::get<0>(deps[d]);
                            int sub_n_segment   = std::get<1>(deps[d]);
                            int sub_n_repeated  = std::get<2>(deps[d]);

                            // Check sub key
                            #ifndef NDEBUG
                            if (d_propagator.find(sub_dep) == d_propagator.end())
                                throw_with_line_number("Could not find sub key '" + sub_dep + "'. ");
                            if (!propagator_finished[sub_dep][sub_n_segment-1])
                                throw_with_line_number("Could not compute '" + key +  "', since '"+ sub_dep + std::to_string(sub_n_segment) + "' is not prepared.");
                            #endif

                            lin_comb<<<N_BLOCKS, N_THREADS>>>(
                                _d_propagator[0], 1.0, _d_propagator[0],
                                sub_n_repeated, d_propagator[sub_dep][sub_n_segment-1], M);
                        }
                        propagator_solver->advance_propagator_discrete(
                            gpu, STREAM,
                            _d_propagator[0],
                            _d_propagator[0],
                            monomer_type, d_q_mask[gpu]);

                        #ifndef NDEBUG
                        propagator_finished[key][0] = true;
                        #endif
                    }
                    else
                    {
                        // Example (four branches)
                        //     A
                        //     |
                        // O - . - B
                        //     |
                        //     C

                        // Legend)
                        // .       : junction
                        // O       : full segment
                        // -, |    : half bonds
                        // A, B, C : other full segments

                        // Combine branches
                        gpu_error_check(cudaMemcpy(d_q_junction[STREAM], d_q_unity[0], sizeof(double)*M, cudaMemcpyDeviceToDevice));

                        for(size_t d=0; d<deps.size(); d++)
                        {
                            std::string sub_dep = std::get<0>(deps[d]);
                            int sub_n_segment   = std::get<1>(deps[d]);

                            // Check sub key
                            #ifndef NDEBUG
                            if (d_propagator.find(sub_dep) == d_propagator.end())
                                throw_with_line_number("Could not find sub key '" + sub_dep + "'. ");
                            if (!propagator_finished[sub_dep][sub_n_segment-1])
                                throw_with_line_number("Could not compute '" + key +  "', since '"+ sub_dep + std::to_string(sub_n_segment) + "' is not prepared.");
                            #endif

                            propagator_solver->advance_propagator_discrete_half_bond_step(
                                gpu, STREAM,
                                d_propagator[sub_dep][sub_n_segment-1],
                                d_q_half_step[STREAM], propagator_analyzer->get_computation_propagator_code(sub_dep).monomer_type);

                            multi_real<<<N_BLOCKS, N_THREADS>>>(d_q_junction[STREAM], d_q_junction[STREAM], d_q_half_step[STREAM], 1.0, M);
                        }
                        gpu_error_check(cudaMemcpy(d_propagator_junction[key], d_q_junction[STREAM], sizeof(double)*M, cudaMemcpyDeviceToDevice));

                        // Add half bond
                        propagator_solver->advance_propagator_discrete_half_bond_step(
                            gpu, STREAM,
                            d_q_junction[STREAM], _d_propagator[0], monomer_type);

                        // Add full segment
                        multi_real<<<N_BLOCKS, N_THREADS>>>(_d_propagator[0], _d_propagator[0], _d_exp_dw, 1.0, M);

                        #ifndef NDEBUG
                        propagator_finished[key][0] = true;
                        #endif
                    }
                }
                else
                {
                    n_segment_from--;
                    is_initialized = false;
                }

                // Multiply mask
                if (is_initialized && d_q_mask[0] != nullptr)
                    multi_real<<<N_BLOCKS, N_THREADS>>>(_d_propagator[0], _d_propagator[0], d_q_mask[0], 1.0, M);

                if (gpu == 0)
                {
                    for(int n=n_segment_from; n<n_segment_to; n++)
                    {
                        #ifndef NDEBUG
                        if (!propagator_finished[key][n-1])
                            throw_with_line_number("unfinished, key: " + key + ", " + std::to_string(n-1));
                        #endif

                        propagator_solver->advance_propagator_discrete(
                            gpu, STREAM, 
                            _d_propagator[n-1],
                            _d_propagator[n],
                            monomer_type, d_q_mask[gpu]);

                        #ifndef NDEBUG
                        propagator_finished[key][n] = true;
                        #endif
                    }
                }
                else if(gpu >= 1)
                {
                    int prev, next;
                    prev = 0;
                    next = 1;

                    // Create events
                    cudaEvent_t kernel_done;
                    cudaEvent_t memcpy_done;
                    gpu_error_check(cudaEventCreate(&kernel_done));
                    gpu_error_check(cudaEventCreate(&memcpy_done));

                    // Copy propagator copy memory from device 1 to device
                    gpu_error_check(cudaMemcpyAsync(
                        d_propagator_device[gpu][prev],
                        _d_propagator[n_segment_from-1],
                        sizeof(double)*M, cudaMemcpyDeviceToDevice, streams[STREAM][1]));

                    gpu_error_check(cudaEventRecord(memcpy_done, streams[STREAM][1]));
                    gpu_error_check(cudaStreamWaitEvent(streams[STREAM][0], memcpy_done, 0));

                    for(int n=n_segment_from; n<n_segment_to; n++)
                    {
                        #ifndef NDEBUG
                        if (!propagator_finished[key][n-1])
                            throw_with_line_number("unfinished, key: " + key + ", " + std::to_string(n-1));
                        #endif

                        // DEVICE 1, STREAM 0: calculate propagators
                        propagator_solver->advance_propagator_discrete(
                            gpu, STREAM, 
                            d_propagator_device[gpu][prev],
                            d_propagator_device[gpu][next],
                            monomer_type, d_q_mask[gpu]);
                        gpu_error_check(cudaEventRecord(kernel_done, streams[STREAM][0]));

                        // DEVICE 1, STREAM 1: copy memory from device 1 to device 0
                        if (n > n_segment_from)
                        {
                            gpu_error_check(cudaMemcpyAsync(
                                _d_propagator[n-1],
                                d_propagator_device[gpu][prev],
                                sizeof(double)*M, cudaMemcpyDeviceToDevice, streams[STREAM][1]));
                            gpu_error_check(cudaEventRecord(memcpy_done, streams[STREAM][1]));
                        }

                        // Wait until computation and memory copy are done
                        gpu_error_check(cudaStreamWaitEvent(streams[STREAM][1], kernel_done, 0));
                        gpu_error_check(cudaStreamWaitEvent(streams[STREAM][0], memcpy_done, 0));

                        std::swap(prev, next);

                        #ifndef NDEBUG
                        propagator_finished[key][n] = true;
                        #endif
                    }

                    // Copy memory from device 1 to device 0
                    gpu_error_check(cudaMemcpyAsync(
                        _d_propagator[n_segment_to-1],
                        d_propagator_device[gpu][prev],
                        sizeof(double)*M, cudaMemcpyDeviceToDevice, streams[STREAM][1]));

                    gpu_error_check(cudaEventRecord(memcpy_done, streams[STREAM][1]));
                    gpu_error_check(cudaStreamWaitEvent(streams[STREAM][0], memcpy_done, 0));
                
                    gpu_error_check(cudaEventDestroy(kernel_done));
                    gpu_error_check(cudaEventDestroy(memcpy_done));
                }
                gpu_error_check(cudaStreamSynchronize(streams[STREAM][0]));
                gpu_error_check(cudaStreamSynchronize(streams[STREAM][1]));
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
            int p                    = std::get<0>(segment_info);
            double *d_propagator_v   = std::get<1>(segment_info);
            double *d_propagator_u   = std::get<2>(segment_info);
            std::string monomer_type = std::get<3>(segment_info);
            int n_aggregated         = std::get<4>(segment_info);
            double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

            single_polymer_partitions[p] = cb->inner_product_inverse_weight_device(
                d_propagator_v,  // q
                d_propagator_u,  // q^dagger
                _d_exp_dw)/n_aggregated/cb->get_volume();
        }

        // Calculate segment concentrations
        for(const auto& block: d_phi_block)
        {
            const auto& key = block.first;
            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            int n_repeated;
            int n_segment_allocated = propagator_analyzer->get_computation_block(key).n_segment_allocated;
            int n_segment_original  = propagator_analyzer->get_computation_block(key).n_segment_original;
            std::string monomer_type = propagator_analyzer->get_computation_block(key).monomer_type;
            double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

            // Contains no '['
            if (dep_u.find('[') == std::string::npos)
                n_repeated = propagator_analyzer->get_computation_block(key).v_u.size();
            else
                n_repeated = 1;

            // Check keys
            #ifndef NDEBUG
            if (d_propagator.find(dep_v) == d_propagator.end())
                throw_with_line_number("Could not find dep_v key'" + dep_v + "'. ");
            if (d_propagator.find(dep_u) == d_propagator.end())
                throw_with_line_number("Could not find dep_u key'" + dep_u + "'. ");
            #endif

            // Calculate phi of one block (possibly multiple blocks when using aggregation)
            calculate_phi_one_block(
                block.second,              // phi
                d_propagator[dep_v],       // dependency v
                d_propagator[dep_u],       // dependency u
                _d_exp_dw,                 // exp_dw
                n_segment_allocated,
                n_segment_original);
            
            // Normalize concentration
            Polymer& pc = molecules->get_polymer(p);
            double norm = molecules->get_ds()*pc.get_volume_fraction()/pc.get_alpha()/single_polymer_partitions[p]*n_repeated;
            lin_comb<<<N_BLOCKS, N_THREADS>>>(block.second, norm, block.second, 0.0, block.second, M);
        }

        // Calculate partition functions and concentrations of solvents
        for(size_t s=0; s<molecules->get_n_solvent_types(); s++)
        {
            double *d_phi_ = d_phi_solvent[s];
            double volume_fraction = std::get<0>(molecules->get_solvent(s));
            std::string monomer_type = std::get<1>(molecules->get_solvent(s));
            double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

            single_solvent_partitions[s] = cb->integral_device(_d_exp_dw)/cb->get_volume();
            linear_scaling_real<<<N_BLOCKS, N_THREADS>>>(d_phi_, _d_exp_dw, volume_fraction/single_solvent_partitions[s], 0.0, M);
        }
        gpu_error_check(cudaSetDevice(0));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaComputationDiscrete::calculate_phi_one_block(
    double *d_phi, double **d_q_1, double **d_q_2, double *d_exp_dw, const int N, const int N_ORIGINAL)
{
    try
    {
        gpu_error_check(cudaSetDevice(0));

        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();

        const int M = cb->get_n_grid();
        // Compute segment concentration
        multi_real<<<N_BLOCKS, N_THREADS>>>(d_phi,d_q_1[N_ORIGINAL-1], d_q_2[0], 1.0, M);
        for(int n=1; n<N; n++)
        {
            add_multi_real<<<N_BLOCKS, N_THREADS>>>(d_phi, d_q_1[N_ORIGINAL-n-1], d_q_2[n], 1.0, M);
        }
        divide_real<<<N_BLOCKS, N_THREADS>>>(d_phi, d_phi, d_exp_dw, 1.0, M);
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
double CudaComputationDiscrete::get_total_partition(int polymer)
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
void CudaComputationDiscrete::get_total_concentration(std::string monomer_type, double *phi)
{
    try
    {
        gpu_error_check(cudaSetDevice(0));

        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();
        const int M = cb->get_n_grid();

        // Initialize to zero
        gpu_error_check(cudaMemset(d_phi, 0, sizeof(double)*M));

        // For each block
        for(const auto& d_block: d_phi_block)
        {
            const auto& key = d_block.first;
            std::string dep_v = std::get<1>(key);
            int n_segment_allocated = propagator_analyzer->get_computation_block(key).n_segment_allocated;
            if (PropagatorCode::get_monomer_type_from_key(dep_v) == monomer_type && n_segment_allocated != 0)
                lin_comb<<<N_BLOCKS, N_THREADS>>>(d_phi, 1.0, d_phi, 1.0, d_block.second, M);
        }

        // For each solvent
        for(int s=0;s<molecules->get_n_solvent_types();s++)
        {
            if (std::get<1>(molecules->get_solvent(s)) == monomer_type)
                lin_comb<<<N_BLOCKS, N_THREADS>>>(d_phi, 1.0, d_phi, 1.0, d_phi_solvent[s], M);
        }
        gpu_error_check(cudaMemcpy(phi, d_phi, sizeof(double)*M, cudaMemcpyDeviceToHost));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaComputationDiscrete::get_total_concentration(int p, std::string monomer_type, double *phi)
{
    try
    {
        gpu_error_check(cudaSetDevice(0));

        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();

        const int M = cb->get_n_grid();
        const int P = molecules->get_n_polymer_types();

        if (p < 0 || p > P-1)
            throw_with_line_number("Index (" + std::to_string(p) + ") must be in range [0, " + std::to_string(P-1) + "]");

        // Initialize to zero
        gpu_error_check(cudaMemset(d_phi, 0, sizeof(double)*M));

        // For each block
        for(const auto& d_block: d_phi_block)
        {
            const auto& key = d_block.first;
            int polymer_idx = std::get<0>(key);
            std::string dep_v = std::get<1>(key);
            int n_segment_allocated = propagator_analyzer->get_computation_block(key).n_segment_allocated;
            if (polymer_idx == p && PropagatorCode::get_monomer_type_from_key(dep_v) == monomer_type && n_segment_allocated != 0)
                lin_comb<<<N_BLOCKS, N_THREADS>>>(d_phi, 1.0, d_phi, 1.0, d_block.second, M);
        }
        gpu_error_check(cudaMemcpy(phi, d_phi, sizeof(double)*M, cudaMemcpyDeviceToHost));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaComputationDiscrete::get_block_concentration(int p, double *phi)
{
    try
    {
        gpu_error_check(cudaSetDevice(0));

        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();

        const int M = cb->get_n_grid();
        const int P = molecules->get_n_polymer_types();

        if (p < 0 || p > P-1)
            throw_with_line_number("Index (" + std::to_string(p) + ") must be in range [0, " + std::to_string(P-1) + "]");

        if (propagator_analyzer->is_aggregated())
            throw_with_line_number("Disable 'aggregation' option to obtain concentration of each block.");

        // Initialize to zero
        gpu_error_check(cudaMemset(d_phi, 0, sizeof(double)*M));

        Polymer& pc = molecules->get_polymer(p);
        std::vector<Block>& blocks = pc.get_blocks();

        for(size_t b=0; b<blocks.size(); b++)
        {
            std::string dep_v = pc.get_propagator_key(blocks[b].v, blocks[b].u);
            std::string dep_u = pc.get_propagator_key(blocks[b].u, blocks[b].v);
            if (dep_v < dep_u)
                dep_v.swap(dep_u);

            lin_comb<<<N_BLOCKS, N_THREADS>>>(d_phi, 0.0, d_phi, 1.0, d_phi_block[std::make_tuple(p, dep_v, dep_u)], M);
            gpu_error_check(cudaMemcpy(&phi[b*M], d_phi, sizeof(double)*M, cudaMemcpyDeviceToHost));
        }
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
double CudaComputationDiscrete::get_solvent_partition(int s)
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
void CudaComputationDiscrete::get_solvent_concentration(int s, double *phi)
{
    try
    {
        gpu_error_check(cudaSetDevice(0));
        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();

        const int M = cb->get_n_grid();
        const int S = molecules->get_n_solvent_types();

        if (s < 0 || s > S-1)
            throw_with_line_number("Index (" + std::to_string(s) + ") must be in range [0, " + std::to_string(S-1) + "]");

        gpu_error_check(cudaMemcpy(phi, d_phi_solvent[s], sizeof(double)*M, cudaMemcpyDeviceToHost));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
std::vector<double> CudaComputationDiscrete::compute_stress()
{
    // This method should be invoked after invoking compute_statistics().

    // To calculate stress, we multiply weighted fourier basis to q(k)*q^dagger(-k).
    // We only need the real part of stress calculation.
    try
    {
        const int N_BLOCKS  = CudaCommon::get_instance().get_n_blocks();
        const int N_THREADS = CudaCommon::get_instance().get_n_threads();
        const int N_GPUS = CudaCommon::get_instance().get_n_gpus();

        const int DIM = cb->get_dim();
        const int M   = cb->get_n_grid();

        std::vector<double> stress(DIM);
        std::map<std::tuple<int, std::string, std::string>, std::array<double,3>> block_dq_dl[N_STREAMS];

        // Reset stress map
        for(const auto& item: d_phi_block)
        {
            for(int i=0; i<N_STREAMS; i++)
                for(int d=0; d<3; d++)
                    block_dq_dl[i][item.first][d] = 0.0;
        }

        // Compute stress for each block
        #pragma omp parallel for num_threads(N_STREAMS) 
        for(size_t b=0; b<d_phi_block.size();b++)
        {
            int gpu;
            const int STREAM = omp_get_thread_num();
            if (N_GPUS == 1)
            {
                gpu = 0;
            }
            else
            {
                gpu = omp_get_thread_num();
            }
            gpu_error_check(cudaSetDevice(gpu));

            auto block = d_phi_block.begin();
            advance(block, b);
            const auto& key   = block->first;

            int p                = std::get<0>(key);
            std::string dep_v    = std::get<1>(key);
            std::string dep_u    = std::get<2>(key);

            const int N           = propagator_analyzer->get_computation_block(key).n_segment_allocated;
            const int N_ORIGINAL  = propagator_analyzer->get_computation_block(key).n_segment_original;
            std::string monomer_type = propagator_analyzer->get_computation_block(key).monomer_type;

            // Contains no '['
            int n_repeated;
            if (dep_u.find('[') == std::string::npos)
                n_repeated = propagator_analyzer->get_computation_block(key).v_u.size();
            else
                n_repeated = 1;

            double **d_q_1 = d_propagator[dep_v];     // Propagator q
            double **d_q_2 = d_propagator[dep_u];     // Propagator q^dagger

            std::array<double,3> _block_dq_dl = {0.0, 0.0, 0.0};

            // Check block_stress_computation_plan
            const auto& _block_stress_compuation_key = block_stress_computation_plan[key];
            if(_block_stress_compuation_key.size() != (unsigned int) (N+1))
            {
                throw_with_line_number("Mismatch of block_stress_computation_plan("
                    + std::to_string(p) + "," + dep_v + "," + dep_u + ") "
                    + std::to_string(_block_stress_compuation_key.size()) + ".size() with N+1 (" + std::to_string(N+1) + ")");
            }

            // Variables for block_stress_computation_plan
            double *d_propagator_v;
            double *d_propagator_u;

            double *d_segment_stress;
            double segment_stress[DIM];
            gpu_error_check(cudaMalloc((void**)&d_segment_stress, sizeof(double)*3));

            int prev, next;
            prev = 0;
            next = 1;

            // Create events
            cudaEvent_t kernel_done;
            cudaEvent_t memcpy_done;
            gpu_error_check(cudaEventCreate(&kernel_done));
            gpu_error_check(cudaEventCreate(&memcpy_done));

            // Copy memory from device to device
            d_propagator_v = std::get<0>(_block_stress_compuation_key[0]);
            d_propagator_u = std::get<1>(_block_stress_compuation_key[0]);

            if (d_propagator_v != nullptr)
            {
                gpu_error_check(cudaMemcpyAsync(&d_q_pair[STREAM][prev][0], d_propagator_v, sizeof(double)*M, cudaMemcpyDeviceToDevice, streams[STREAM][1]));
                gpu_error_check(cudaMemcpyAsync(&d_q_pair[STREAM][prev][M], d_propagator_u, sizeof(double)*M, cudaMemcpyDeviceToDevice, streams[STREAM][1]));
                gpu_error_check(cudaEventRecord(memcpy_done, streams[STREAM][1]));
            }
            gpu_error_check(cudaStreamWaitEvent(streams[STREAM][0], memcpy_done, 0));

            for(int n=0; n<=N; n++)
            {
                // STREAM 1: copy memory from device to device
                if (n+1 <= N)
                {
                    d_propagator_v = std::get<0>(_block_stress_compuation_key[n+1]);
                    d_propagator_u = std::get<1>(_block_stress_compuation_key[n+1]);

                    if (d_propagator_v != nullptr)
                    {
                        gpu_error_check(cudaMemcpyAsync(&d_q_pair[STREAM][next][0], d_propagator_v, sizeof(double)*M, cudaMemcpyDeviceToDevice, streams[STREAM][1]));
                        gpu_error_check(cudaMemcpyAsync(&d_q_pair[STREAM][next][M], d_propagator_u, sizeof(double)*M, cudaMemcpyDeviceToDevice, streams[STREAM][1]));
                        gpu_error_check(cudaEventRecord(memcpy_done, streams[STREAM][1]));
                    }
                }

                // STREAM 0: Compute stress
                d_propagator_v = std::get<0>(_block_stress_compuation_key[n]);
                bool is_half_bond_length = std::get<2>(_block_stress_compuation_key[n]);
                if (d_propagator_v != nullptr)
                {
                    propagator_solver->compute_single_segment_stress_discrete(
                        gpu, STREAM, d_q_pair[STREAM][prev], d_segment_stress, monomer_type, is_half_bond_length);
                    gpu_error_check(cudaEventRecord(kernel_done, streams[STREAM][0]));
                }

                // Wait until computation and memory copy are done
                gpu_error_check(cudaStreamWaitEvent(streams[STREAM][1], kernel_done, 0));
                gpu_error_check(cudaStreamWaitEvent(streams[STREAM][0], memcpy_done, 0));

                if (d_propagator_v != nullptr)
                {
                    gpu_error_check(cudaMemcpy(segment_stress, d_segment_stress, sizeof(double)*DIM, cudaMemcpyDeviceToHost));
                    for(int d=0; d<DIM; d++)
                        _block_dq_dl[d] += segment_stress[d]*n_repeated;
                }
                std::swap(prev, next);
            }
            gpu_error_check(cudaStreamSynchronize(streams[STREAM][0]));
            gpu_error_check(cudaStreamSynchronize(streams[STREAM][1]));
            gpu_error_check(cudaEventDestroy(kernel_done));
            gpu_error_check(cudaEventDestroy(memcpy_done));

            // Copy stress data
            for(int d=0; d<DIM; d++)
                block_dq_dl[STREAM][key][d] += _block_dq_dl[d];
                
            cudaFree(d_segment_stress);
        }
        // Synchronize all GPUs
        for(int gpu=0; gpu<N_GPUS; gpu++)
        {
            gpu_error_check(cudaSetDevice(gpu));
            gpu_error_check(cudaDeviceSynchronize());
        }

        gpu_error_check(cudaSetDevice(0));
        // Compute total stress
        for(int d=0; d<DIM; d++)
            stress[d] = 0.0;
        for(const auto& block: d_phi_block)
        {
            const auto& key = block.first;
            int p             = std::get<0>(key);
            std::string dep_v = std::get<1>(key);
            std::string dep_u = std::get<2>(key);
            Polymer& pc  = molecules->get_polymer(p);

            for(int i=0; i<N_STREAMS; i++)
                for(int d=0; d<DIM; d++)
                    stress[d] += block_dq_dl[i][key][d]*pc.get_volume_fraction()/pc.get_alpha()/single_polymer_partitions[p];
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
void CudaComputationDiscrete::get_chain_propagator(double *q_out, int polymer, int v, int u, int n)
{ 
    // This method should be invoked after invoking compute_statistics()

    // Get chain propagator for a selected polymer, block and direction.
    // This is made for debugging and testing.
    try
    {
        const int M = cb->get_n_grid();
        Polymer& pc = molecules->get_polymer(polymer);
        std::string dep = pc.get_propagator_key(v,u);

        if (propagator_analyzer->get_computation_propagator_codes().find(dep) == propagator_analyzer->get_computation_propagator_codes().end())
            throw_with_line_number("Could not find the propagator code '" + dep + "'. Disable 'aggregation' option to obtain propagator_analyzer.");

        const int N = propagator_analyzer->get_computation_propagator_codes()[dep].max_n_segment;
        if (n < 1 || n > N)
            throw_with_line_number("n (" + std::to_string(n) + ") must be in range [1, " + std::to_string(N) + "]");

        gpu_error_check(cudaMemcpy(q_out, d_propagator[dep][n-1], sizeof(double)*M,cudaMemcpyDeviceToHost));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
bool CudaComputationDiscrete::check_total_partition()
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
    for(const auto& block: d_phi_block)
    {
        const auto& key = block.first;
        int p                = std::get<0>(key);
        std::string dep_v    = std::get<1>(key);
        std::string dep_u    = std::get<2>(key);

        int n_aggregated;
        int n_segment_allocated = propagator_analyzer->get_computation_block(key).n_segment_allocated;
        int n_segment_original  = propagator_analyzer->get_computation_block(key).n_segment_original;
        std::string monomer_type = propagator_analyzer->get_computation_block(key).monomer_type;
        double *_d_exp_dw = propagator_solver->d_exp_dw[0][monomer_type];

        // std::cout<< p << ", " << dep_v << ", " << dep_u << ": " << n_segment_original << ", " << n_segment_allocated << std::endl;

        // Contains no '['
        if (dep_u.find('[') == std::string::npos)
            n_aggregated = 1;
        else
            n_aggregated = propagator_analyzer->get_computation_block(key).v_u.size();

        for(int n=0;n<n_segment_allocated;n++)
        {
            double total_partition = cb->inner_product_inverse_weight_device(
                d_propagator[dep_v][(n_segment_original-n-1)],
                d_propagator[dep_u][n], _d_exp_dw)/n_aggregated/cb->get_volume();

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