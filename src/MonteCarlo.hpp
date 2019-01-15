#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"


#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

class MonteCarlo
{
public:
    BlackScholesModel *mod_; /*! pointeur vers le modèle */
    Option *opt_; /*! pointeur sur l'option */
    PnlRng *rng_; /*! pointeur sur le générateur */
    double fdStep_; /*! pas de différence finie */
    int nbSamples_; /*! nombre de tirages Monte Carlo */

    MonteCarlo();

    MonteCarlo(BlackScholesModel *mod, Option *opt, PnlRng *rng, double fdStep, int nbSamples);

    //MonteCarlo(const MonteCarlo & mc);

    virtual ~MonteCarlo();

    /**
     * Calcule le prix de l'option à la date 0, version maître
     *
     * @param[out] prix valeur de l'estimateur Monte Carlo
     * @param[out] ic largeur de l'intervalle de confiance
     */
    void price_master(double &prix, double &ic);

    /**
     * Calcule le prix de l'option à la date 0, version maître,
     * avec largeur d'intervalle de confiance recherchée fixée
     *
     * @param[out] prix valeur de l'estimateur Monte Carlo
     * @param[out] ic largeur de l'intervalle de confiance
     */
    void price_master_given_precision(double &prix, double &ic, double precision);

    /**
     * Calcule le prix de l'option à la date 0, version esclave
     *
     */
    void price_slave();

    /**
     * Calcule le prix de l'option à la date 0, version esclave,
     * avec largeur d'intervalle de confiance recherchée fixée
     *
     */
    void price_slave_given_precision();
};
