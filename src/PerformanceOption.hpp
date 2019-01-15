
#ifndef PERFORMANCEOPTION_H
#define PERFORMANCEOPTION_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

#include "Option.hpp"

class PerformanceOption : public Option
{
private:
    PnlVect* coefs_;
    PerformanceOption();
    PnlVect* payoffVectMemSpaceSums_; //sums
    PnlVect* payoffVectMemSpaceNumerateur_; //numerateur
    PnlVect* payoffVectMemSpaceDenominateur_; //denominateur
public:
    PerformanceOption(double T, double nbTimeSteps, double size, PnlVect* coefs, double strike);
    ~PerformanceOption();
    double payoff(const PnlMat *path);
};
#endif
