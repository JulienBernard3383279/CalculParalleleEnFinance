#ifndef BASKETOPTION_H
#define BASKETOPTION_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

#include "Option.hpp"

class BasketOption : public Option
{
private:
    PnlVect* coefs_;
    double strike_;
    PnlVect* payoffVectMemSpace_;
    BasketOption();

public:
    BasketOption(double T, double nbTimeSteps, double size, PnlVect* coefs, double strike);
    ~BasketOption();

    double payoff(const PnlMat *path);
};

#endif
