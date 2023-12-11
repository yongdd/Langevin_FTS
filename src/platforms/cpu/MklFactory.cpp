/*----------------------------------------------------------
* class MklFactory
*-----------------------------------------------------------*/

#include <iostream>
#include <array>
#include <vector>
#include <string>

#include "mkl.h"

#include "CpuArray.h"
#include "CpuComputationBox.h"
#include "CpuSolverContinuous.h"
#include "CpuSolverDiscrete.h"
#include "CpuAndersonMixing.h"
#include "MklFactory.h"

MklFactory::MklFactory(bool reduce_memory_usage)
{
    this->reduce_memory_usage = reduce_memory_usage;

    if (this->reduce_memory_usage)
        std::cout << "(warning) Reducing memory usage option only works for CUDA. This option will be ignored in MKL." << std::endl;

}
Array* MklFactory::create_array(
    unsigned int size)
{
    return new CpuArray(size);
}

Array* MklFactory::create_array(
    double *data,
    unsigned int size)
{
    return new CpuArray(data, size);
}
ComputationBox* MklFactory::create_computation_box(
    std::vector<int> nx, std::vector<double> lx, const double *mask)
{
    return new CpuComputationBox(nx, lx, mask);
}
Molecules* MklFactory::create_molecules_information(
    std::string chain_model, double ds, std::map<std::string, double> bond_lengths) 
{
    return new Molecules(chain_model, ds, bond_lengths);
}
Solver* MklFactory::create_pseudospectral_solver(ComputationBox *cb, Molecules *molecules, PropagatorsAnalyzer* propagators_analyzer)
{
    std::string chain_model = molecules->get_model_name();
    if ( chain_model == "continuous" )
    {
        return new CpuSolverContinuous(cb, molecules, propagators_analyzer);
    }
    else if ( chain_model == "discrete" )
    {
        return new CpuSolverDiscrete(cb, molecules, propagators_analyzer);
    }
    return NULL;
}
AndersonMixing* MklFactory::create_anderson_mixing(
    int n_var, int max_hist, double start_error,
    double mix_min, double mix_init)
{
    return new CpuAndersonMixing(
        n_var, max_hist, start_error, mix_min, mix_init);
}
void MklFactory::display_info()
{
    MKLVersion Version;
 
    mkl_get_version(&Version);
    std::cout<< "==================== MKL Version ====================" << std::endl;
    printf("Major version:           %d\n",Version.MajorVersion);
    printf("Minor version:           %d\n",Version.MinorVersion);
    printf("Update version:          %d\n",Version.UpdateVersion);
    printf("Product status:          %s\n",Version.ProductStatus);
    printf("Build:                   %s\n",Version.Build);
    printf("Platform:                %s\n",Version.Platform);
    printf("Processor optimization:  %s\n",Version.Processor);
    printf("================================================================\n");
}
