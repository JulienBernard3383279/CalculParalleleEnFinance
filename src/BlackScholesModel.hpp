#ifndef BLACKSCHOLESMODEL_H
#define BLACKSCHOLESMODEL_H

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Modèle de Black Scholes
class BlackScholesModel
{
    PnlMat *gammaMemSpace_;
    PnlMat *gMemSpace_;
    PnlVect tempMemSpace1_;
    PnlVect tempMemSpace2_;
public:
    int size_; /// nombre d'actifs du modèle
    double r_; /// taux d'intérêt
    double rho_; /// paramètre de corrélation
    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *spot_; /// valeurs initiales du sous-jacent

    PnlVect *trend_; /// tendance du modèle

    BlackScholesModel();
    BlackScholesModel(int size, double r, double rho, PnlVect *sigma, 
                      PnlVect *spot, PnlVect *trend);
    BlackScholesModel(const BlackScholesModel & bsm);
    virtual ~BlackScholesModel();


    /**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (nbTimeSteps+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     * @param[in] rng générateur de nombres aléatoires
     */
    void asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);

    /**
     * Dans le cas d'appel d'asset successifs avec la configuration actuelle,
     * initialise les zones mémorielles dédiées à asset avec les matrices ou
     * vecteurs de bonnes dimensions.
     * 
     *  @param[in] nbTimeSteps nombre de dates de constatation
     */
    void initAsset(int nbTimeSteps);

    /**
     * Génère une trajectoire du modèle et la stocke dans path. Les appels doivent être fait après un initAsset.
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (nbTimeSteps+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     * @param[in] rng générateur de nombres aléatoires
     */
    void postInitAsset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);

    /**
     * @brief Shift d'une trajectoire du sous-jacent
     *
     * @param[in]  path contient en input la trajectoire
     * du sous-jacent
     * @param[out] shift_path contient la trajectoire path
     * dont la composante d a été shiftée par (1+h)
     * à partir de la date t.
     * @param[in] t date à partir de laquelle on shift
     * @param[in] h pas de différences finies
     * @param[in] d indice du sous-jacent à shifter
     * @param[in] timestep pas de constatation du sous-jacent
     */
     void shiftAsset(PnlMat *shift_path, const PnlMat *path, int d,
            double h, double t, double timestep);

    /**
     * Génère une simulation du marché
     *
     * @param[out] market contient une trajectoire du modèle.
     * C'est une matrice de taille (H+1) x d
     * @param[in] T maturité
     * @param[in] H nombre de dates de constatation
     * @param[in] rng générateur de nombres aléatoires
     */
    void simul_market(PnlMat *market, double T, int H, PnlRng *rng);
};

#endif
