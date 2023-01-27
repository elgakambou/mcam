#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

/// \brief Performance Option
class PerformanceOption : public Option
{
public:
    double strike_; // strike > 0
    PnlVect* lambda_; // lambda, un vecteur de réels
    PerformanceOption(double T, int dates, int size, double strike, PnlVect* lambda);
    double payoff(const PnlMat *path, double t);
};
