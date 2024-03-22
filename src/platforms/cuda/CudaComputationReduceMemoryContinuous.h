/*-------------------------------------------------------------
This is a derived CudaComputationReduceMemoryContinuous class

GPU memory usage is reduced by storing propagators in main memory.
In the GPU memory, array space that can store only two steps of propagator is allocated.
There are two streams. One is responsible for data transfers between CPU and GPU, another is responsible
for kernel executions. Overlapping of kernel execution and data transfers is utilized so that 
they can be simultaneously executed. As a result, data transfer time can be hided. For more explanation,
please see the supporting information of [Macromolecules 2021, 54, 24, 11304].
*------------------------------------------------------------*/

#ifndef CUDA_PSEUDO_REDUCE_MEMORY_CONTINUOUS_H_
#define CUDA_PSEUDO_REDUCE_MEMORY_CONTINUOUS_H_

#include <array>
#include <cufft.h>

#include "ComputationBox.h"
#include "Polymer.h"
#include "Molecules.h"
#include "PropagatorComputation.h"
#include "CudaCommon.h"
#include "CudaSolver.h"
#include "Scheduler.h"

class CudaComputationReduceMemoryContinuous : public PropagatorComputation
{
private:
    // Pseudo-spectral PDE solver
    CudaSolver *propagator_solver;

    std::string method;

    // The number of parallel streams
    static const int N_STREAMS = 2;

    // Two streams for each gpu
    cudaStream_t streams[N_STREAMS][2]; // one for kernel execution, the other for memcpy

    // For pseudo-spectral: advance_one propagator()
    double *d_q_one[N_STREAMS][2];               // one for prev, the other for next
    double *d_propagator_sub_dep[N_STREAMS][2];  // one for prev, the other for next

    // All elements are 1 for initializing propagators
    double *d_q_unity[MAX_GPUS];

    // q_mask to make impenetrable region for nano particles
    double *d_q_mask[MAX_GPUS];

    // For concentration computation
    double *d_q_block_v[2];    // one for prev, the other for next
    double *d_q_block_u[2];    // one for prev, the other for next
    double *d_phi;

    // For stress calculation: compute_stress()
    double *d_q_pair[N_STREAMS][2];  // one for prev, the other for next

    // Scheduler for propagator computation 
    Scheduler *sc;

    // Host pinned memory space to store propagator, key: (dep) + monomer_type, value: propagator
    std::map<std::string, double **> propagator;
    // Map for deallocation of d_propagator
    std::map<std::string, int> propagator_size;
    // Check if computation of propagator is finished
    #ifndef NDEBUG
    std::map<std::string, bool *> propagator_finished;
    #endif

    // Total partition function
    double *single_polymer_partitions; 
    // Remember one segment for each polymer chain to compute total partition function
    // (polymer id, propagator forward, propagator backward, n_repeated)
    std::vector<std::tuple<int, double *, double *, int>> single_partition_segment;

    // Host pinned space to store concentration, key: (polymer id, dep_v, dep_u) (assert(dep_v <= dep_u)), value: concentration
    std::map<std::tuple<int, std::string, std::string>, double *> phi_block;

    // Total partition functions for each solvent
    double* single_solvent_partitions;

    // Solvent concentrations
    std::vector<double *> phi_solvent;

    // Calculate concentration of one block
    void calculate_phi_one_block(double *phi, double **q_1, double **q_2, const int N, const int N_OFFSET, const double NORM);

    // Compute statistics with inputs from selected device arrays
    void compute_statistics(std::string device,
        std::map<std::string, const double*> w_input,
        std::map<std::string, const double*> q_init = {});
public:

    CudaComputationReduceMemoryContinuous(ComputationBox *cb, Molecules *pc, PropagatorAnalyzer *propagator_analyzer, std::string method);
    ~CudaComputationReduceMemoryContinuous();

    void update_laplacian_operator() override;
    void compute_statistics(
        std::map<std::string, const double*> w_block,
        std::map<std::string, const double*> q_init = {}) override
    {
        compute_statistics("cpu", w_block, q_init);
    };
    void compute_statistics_device(
        std::map<std::string, const double*> d_w_block,
        std::map<std::string, const double*> d_q_init = {}) override
    {
        compute_statistics("gpu", d_w_block, d_q_init);
    };
    double get_total_partition(int polymer) override;
    void get_total_concentration(std::string monomer_type, double *phi) override;
    void get_total_concentration(int polymer, std::string monomer_type, double *phi) override;
    void get_block_concentration(int polymer, double *phi) override;
    std::vector<double> compute_stress() override;
    void get_chain_propagator(double *q_out, int polymer, int v, int u, int n) override;

    double get_solvent_partition(int s) override;
    void get_solvent_concentration(int s, double *phi) override;

    // For tests
    bool check_total_partition() override;
};

#endif
