#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"

/// \brief Basket Option
class BasketOption : public Option
{
public:
    double strike_; // strike, un réel (le cas K <= correspond à un put)
    PnlVect* lambda_; // lambda, un vecteur de réels
    BasketOption(double T, int dates, int size, double strike, PnlVect* lambda);
    double payoff(const PnlMat *path, double t);
};


