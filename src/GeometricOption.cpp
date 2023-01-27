#include "GeometricOption.hpp"
#include <cmath>

using namespace std;

GeometricOption::GeometricOption(double T, int dates, int size, double strike)
{
    T_ = T;
    dates_ = dates;
    size_ = size;
    strike_ = strike;
}

double GeometricOption::payoff(const PnlMat *path, double t)
{
    int tk = (int)(t * dates_ / T_);
    PnlVect *Stk = pnl_vect_create(size_);
    pnl_mat_get_col(Stk, path, tk);
    double prod = pnl_vect_prod(Stk);
    pnl_vect_free(&Stk);
    return max(strike_ - pow(prod, (double)1.0 / ((double)size_)), 0.0);
}
