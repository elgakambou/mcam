#pragma once

#include "pnl/pnl_basis.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

class Quantifier {
    public :
        int id;
        int sizeOfX;
        PnlBasis* B;
        virtual void estimate ( PnlMat* dataX, PnlVect* dataY, PnlVect* coeff) = 0;
        virtual  double eval(const PnlVect* coeff, const PnlVect* x) = 0;
};