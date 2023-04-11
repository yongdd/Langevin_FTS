/*-------------------------------------------------------------
* This is an abstract AndersonMixing class
*------------------------------------------------------------*/

#ifndef ANDERSON_MIXING_H_
#define ANDERSON_MIXING_H_

#include <cassert>
#include <iostream>

#include "Exception.h"

class AndersonMixing
{
protected:
    int n_var, max_hist, n_anderson;
    double start_error, mix_min, mix, mix_init;

    void find_an(double **u, double *v, double *a, int n);
public:
    AndersonMixing(int n_var, int max_hist, double start_error, double mix_min, double mix_init);
    virtual ~AndersonMixing(){};

    virtual void reset_count(){};
    int get_n_var(){ return n_var;};
    virtual void calculate_new_fields(
        double *w_new, double *w_current, double *w_deriv,
        double old_error_level, double error_level)=0;
};
#endif
