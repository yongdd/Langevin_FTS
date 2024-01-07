/*-------------------------------------------------------------
* This is a derived CudaComputationContinuous class
*------------------------------------------------------------*/

#ifndef CUDA_COMPUTATION_CONTINUOUS_H_
#define CUDA_COMPUTATION_CONTINUOUS_H_

#include <array>
#include <cufft.h>

#include "ComputationBox.h"
#include "Polymer.h"
#include "Molecules.h"
#include "PropagatorComputation.h"
#include "CudaCommon.h"
#include "CudaSolver.h"
#include "Scheduler.h"

class CudaComputationContinuous : public PropagatorComputation
{
private:
    // Pseudo-spectral PDE solver
    CudaSolver *propagator_solver;

    std::string method;

    // The number of parallel streams
    static const int N_STREAMS = 2;

    // Two streams for each gpu
    cudaStream_t streams[N_STREAMS][2]; // one for kernel execution, the other for memcpy

    // All elements are 1 for initializing propagators
    double *d_q_unity[MAX_GPUS]; 

    // q_mask to make impenetrable region for nano particles
    double *d_q_mask[MAX_GPUS];

    // One for prev, the other for next
    double *d_q_pair[N_STREAMS][2];

    // Scheduler for propagator computation 
    Scheduler *sc;
    // gpu memory space to store propagator, key: (dep) + monomer_type, value: propagator
    std::map<std::string, double **> d_propagator; 
    // Map for deallocation of d_propagator
    std::map<std::string, int> propagator_size;
    // Temporary arrays for device_1, one for prev, the other for next
    double *d_propagator_device[MAX_GPUS][2];
    // Check if computation of propagator is finished
    #ifndef NDEBUG  
    std::map<std::string, bool *> propagator_finished;
    #endif

    // Total partition functions
    double *single_polymer_partitions; 
    // Remember one segment for each polymer chain to compute total partition function
    // (polymer id, propagator forward, propagator backward, n_aggregated)
    std::vector<std::tuple<int, double *, double *, int>> single_partition_segment;

    // gpu memory space to store concentration, key: (polymer id, dep_v, dep_u) (assert(dep_v <= dep_u)), value: concentration
    std::map<std::tuple<int, std::string, std::string>, double *> d_phi_block;
    // Temp array for concentration computation
    double *d_phi;

    // Total partition functions for each solvent
    double* single_solvent_partitions;

    // Solvent concentrations
    std::vector<double *> d_phi_solvent;

    // Calculate concentration of one block
    void calculate_phi_one_block(double *d_phi, double **d_q_1, double **d_q_2, const int N, const int N_OFFSET, const int N_ORIGINAL);

    // Compute statistics with inputs from selected device arrays
    void compute_statistics(std::string device,
        std::map<std::string, const double*> w_input,
        std::map<std::string, const double*> q_init = {});
public:

    CudaComputationContinuous(ComputationBox *cb, Molecules *pc, PropagatorAnalyzer *propagator_analyzer, std::string method);
    ~CudaComputationContinuous();

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
