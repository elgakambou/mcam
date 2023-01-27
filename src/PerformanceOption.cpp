#include "PerformanceOption.hpp"
#include <cmath>

using namespace std;

PerformanceOption::PerformanceOption(double T, int dates, int size, double strike, PnlVect* lambda)
{
    T_ = T;
    dates_ = dates;
    size_ = size;
    strike_ = strike;
    lambda_ = lambda;
}

double PerformanceOption::payoff(const PnlMat *path, double t)
{
    int tk = (int)(t * dates_ / T_);
    PnlVect *Stk = pnl_vect_create(size_);
    pnl_mat_get_col(Stk, path, tk);
    pnl_vect_mult_vect_term(Stk, lambda_);
    double Stk_max = pnl_vect_max(Stk);
    pnl_vect_free(&Stk);
    return max(Stk_max - strike_, 0.0);
}
