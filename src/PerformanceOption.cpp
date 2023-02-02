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

double PerformanceOption::payoff(const PnlVect *spots)
{
    PnlVect *prod = pnl_vect_copy(spots);
    pnl_vect_mult_vect_term(prod, lambda_);
    double payoff = max(pnl_vect_max(prod) - strike_, 0.0);
    pnl_vect_free(&prod);
    return payoff;
}
