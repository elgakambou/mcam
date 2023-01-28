#pragma once

#include "Quantifier.hpp"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

class PolynomialRegression :  public Quantifier {

    public :
        int maxDegree;
        PolynomialRegression(int maxDegree, int sizeOfX);
        void estimate ( PnlMat* dataX, PnlVect* dataY, PnlVect* coeff);
        double eval(const PnlVect* coeff, const PnlVect* x);
    private :
        PnlMat* m_dataX;
        PnlVect* m_dataY;

};