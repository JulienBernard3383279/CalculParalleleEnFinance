#include <exception>
#include <iostream>
#include <cstring>
#include <string>
#include <time.h>
#include <stdio.h>

#include "mpi.h"

#include "Parser.hpp"
#include "Option.hpp"
#include "MonteCarlo.hpp"
#include "BasketOption.hpp"
#include "AsianOption.hpp"
#include "PerformanceOption.hpp"
#include "BlackScholesModel.hpp"

using namespace std;

int main (int argc, char *argv[]) {
    int mpiSize, mpiRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpiRank);

    /* TEMPS DE CALCUL */
    double calcTime;
    clock_t start, end;
    
    /* BLACK-SCHOLES */
    int size; /// nombre d'actifs du modèle
    double r; /// taux d'intérêt
    double rho; /// paramètre de corrélation
    PnlVect *sigma; /// vecteur de volatilités
    PnlVect *spot; /// valeurs initiales du sous-jacent
    double T;

    /* OPTION */
    double strike;
    PnlVect *coefs;
    int nbTimeSteps;
    double fdStep = 1; /*! pas de différence finie */
    size_t nbSamples; /*! nombre de tirages Monte Carlo */
    string optionType;
    PnlVect *trend; /// tendance du modèle

    // Option de l'exécutable
    string exeOption;

    // Parser
    Param *P;

    /* Lecture des arguments */
    if (argc != 2 && argc != 3) {
        cout << "Utilisation : ./parser data_input [precision]" << endl;
        return EXIT_FAILURE;
    }

    bool precisionMode = (argc == 3);
    double precision=0;
    if (precisionMode) {
        precision=atof(argv[2]);
    }

    string data_input = argv[1];
    P = new Parser(data_input.c_str());

    /* Parsing */
    P->extract("option type", optionType);
    P->extract("maturity", T);
    P->extract("option size", size);
    if (optionType != "performance") { P->extract("strike", strike); }
    P->extract("payoff coefficients", coefs, size);
    P->extract("sample number", nbSamples);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("interest rate", r);
    P->extract("correlation", rho);
    P->extract("timestep number", nbTimeSteps);
    P->extract("trend", trend, size);
    if (trend==NULL) {
        trend = pnl_vect_create_from_zero(size);
    }
    
    /* Caractéristiques de l'option */
    if (!mpiRank) {
        cout << endl << "option size\t\t\t" << size << endl;
        if (optionType != "performance") { cout << "strike\t\t\t\t" << strike << endl; }
        cout << "spot\t\t\t\t";
        pnl_vect_print_asrow(spot);
        cout << "maturity\t\t\t" << T << endl;
        cout << "volatility\t\t\t";
        pnl_vect_print_asrow(sigma);
        cout << "interest rate\t\t\t" << r << endl;
        cout << "correlation\t\t\t" << rho << endl;
        cout << "trend\t\t\t\t";
        pnl_vect_print_asrow(trend);

        cout << endl << "option type\t\t\t" << optionType << endl;
        cout << "payoff coefficients\t\t";
        pnl_vect_print_asrow(coefs);

        cout << endl << "timestep number\t\t\t" << nbTimeSteps << endl;
        cout << "sample number\t\t\t" << nbSamples << endl;
    }
    /* Rng et BS Model en t = 0 */
    PnlRng *rng = pnl_rng_dcmt_create_id(mpiRank, 0);
    pnl_rng_sseed(rng, time(NULL)); //mpiRank = le mal
    BlackScholesModel *mod = new BlackScholesModel(size, r, rho, sigma, spot, trend);

    /* Création de path */
    PnlMat *path;
    
    path = pnl_mat_create(nbTimeSteps+1, size);
    //mod->asset(path, 0, nbTimeSteps, rng);

    // MonteCarlo
    MonteCarlo *mc;

    /* Construction de l'option */
    Option *opt;
    if (optionType == "basket") {
        opt = new BasketOption(T, nbTimeSteps, size, coefs, strike);
        mc = new MonteCarlo(mod, opt, rng, fdStep, nbSamples);
    }
    else if (optionType == "asian") {
        opt = new AsianOption(T, nbTimeSteps, size, coefs, strike);
        mc = new MonteCarlo(mod, opt, rng, fdStep, nbSamples);
    }
    else if (optionType == "performance") {
        opt = new PerformanceOption(T, nbTimeSteps, size, coefs, strike);
        mc = new MonteCarlo(mod, opt, rng, fdStep, nbSamples);
    }
    else {
        cerr << "Erreur sur le type de l'option" << endl;
        return EXIT_FAILURE;
    }

    if (mpiRank) {
        precisionMode ? mc -> price_slave_given_precision() : mc->price_slave();
    }
    else {
        /* Calcul du prix  en t = 0 */
        cout << endl << "***** Valeurs calculées en t = 0 *****" << endl;

        double prix, ic;
        start = clock();
        precisionMode ? mc-> price_master_given_precision(prix, ic, precision) : mc->price_master(prix, ic);
        end = clock();

        cout << "Prix de l'option :\t\t" << prix << endl;
        cout << "Largeur de l'intervalle de confiance à 95pct:\t" << ic << endl;
        cout << "Ecart type :\t\t\t" << ic / 4.0 << endl;
        cout << endl;
        cout << "Temps de calcul : \t\t" << ((float)(end - start))/CLOCKS_PER_SEC << " secondes" << endl;
        cout << endl;
    }

    pnl_rng_free(&rng);
    pnl_vect_free(&spot);
    pnl_vect_free(&sigma);
    pnl_vect_free(&trend);
    pnl_vect_free(&coefs);
    pnl_mat_free(&path);
    
    delete P;
    delete mc;
    delete mod;
    delete opt;

    MPI_Finalize();
    return EXIT_SUCCESS;
}
