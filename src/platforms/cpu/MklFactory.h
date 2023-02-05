/*----------------------------------------------------------
* class MklFactory
*-----------------------------------------------------------*/

#ifndef MKL_FACTORY_H_
#define MKL_FACTORY_H_

#include "ComputationBox.h"
#include "PolymerChain.h"
#include "Mixture.h"
#include "Pseudo.h"
#include "AndersonMixing.h"
#include "AbstractFactory.h"

class MklFactory : public AbstractFactory
{
public :
    MklFactory(std::string chain_model);

    ComputationBox* create_computation_box(
        std::vector<int> nx,
        std::vector<double> lx) override;

    Mixture* create_mixture(
        double ds, std::map<std::string, double> bond_lengths, bool use_superposition) override;

    Pseudo* create_pseudo(ComputationBox *cb, Mixture *mx) override;

    AndersonMixing* create_anderson_mixing(
        int n_var, int max_hist, double start_error,
        double mix_min, double mix_init) override;

    void display_info() override;
};
#endif
