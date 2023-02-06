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

double GeometricOption::payoff(const PnlVect *spots)
{
    return max(strike_ - pow(pnl_vect_prod(spots), 1.0 / ((double)size_)), 0.0);
}
