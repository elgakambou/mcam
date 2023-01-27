#pragma once

#include "Option.hpp"

/// \brief Geometric Option (put)
class GeometricOption : public Option
{
public:
    double strike_; // strike > 0
    GeometricOption(double T, int dates, int size, double strike);
    double payoff(const PnlMat *path, double t);
};


