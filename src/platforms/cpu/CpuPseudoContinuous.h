/*----------------------------------------------------------
* This class defines a derived class for pseudo-spectral method
*-----------------------------------------------------------*/

#ifndef CPU_PSEUDO_CONTINUOUS_H_
#define CPU_PSEUDO_CONTINUOUS_H_

#include <string>
#include <vector>
#include <map>

#include "Exception.h"
#include "ComputationBox.h"
#include "Pseudo.h"

class CpuPseudoContinuous : public Pseudo
{
private:

public:
    CpuPseudoContinuous(ComputationBox *cb);
    ~CpuPseudoContinuous();
    void update_bond_function() override;
};
#endif
