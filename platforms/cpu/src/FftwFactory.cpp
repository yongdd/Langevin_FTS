/*----------------------------------------------------------
* class FftwFactory
*-----------------------------------------------------------*/

#include <array>
#include <vector>
#include <string>
#include <algorithm>

#include "FftwFFT3D.h"
#include "CpuPseudoGaussian.h"
#include "CpuPseudoDiscrete.h"
#include "CpuAndersonMixing.h"
#include "FftwFactory.h"

PolymerChain* FftwFactory::create_polymer_chain(double f, int NN, double chi_n)
{
    return new PolymerChain(f, NN, chi_n);
}
SimulationBox* FftwFactory::create_simulation_box(
    std::array<int,3> nx, std::array<double,3>  lx)
{
    return new SimulationBox(nx, lx);
}
Pseudo* FftwFactory::create_pseudo(SimulationBox *sb, PolymerChain *pc, std::string str_model)
{
    std::transform(str_model.begin(), str_model.end(), str_model.begin(),
    [](unsigned char c){ return std::tolower(c); });
    
    if ( str_model == "gaussian" )
        return new CpuPseudoGaussian(sb, pc, new FftwFFT3D(sb->get_nx()));
    else if ( str_model == "discrete" )
        return new CpuPseudoDiscrete(sb, pc, new FftwFFT3D(sb->get_nx()));
    return NULL;
}
AndersonMixing* FftwFactory::create_anderson_mixing(
    SimulationBox *sb, int n_comp,
    double max_anderson, double start_anderson_error,
    double mix_min, double mix_init)
{
    return new CpuAndersonMixing(
               sb, n_comp, max_anderson,
               start_anderson_error, mix_min, mix_init);
}
