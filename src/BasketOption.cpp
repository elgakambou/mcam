#include "BasketOption.hpp"
#include <cmath>

using namespace std;

BasketOption::BasketOption(double T, int dates, int size, double strike, PnlVect* lambda)
{
    T_ = T;
    dates_ = dates;
    size_ = size;
    strike_ = strike;
    lambda_ = lambda;
}

double BasketOption::payoff(const PnlMat *path, double t)
{
    int tk = (int)(t * dates_ / T_);
    PnlVect *Stk = pnl_vect_create(size_);
    pnl_mat_get_col(Stk, path, tk);
    double sum = pnl_vect_scalar_prod(lambda_, Stk);
    pnl_vect_free(&Stk);
    return max(sum - strike_, 0.0);
}
