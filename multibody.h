#include <armadillo>

// Class that handles perturbation theory of a polymer
class GaussSystem {
public:
    int N;
    arma::mat lap, lti, rg2, opm;
    double rg20, D, a;
 
    double _z0 = -1.0;
 
    GaussSystem(arma::mat laplacian, double dimension, double seglen);
    double term1(); //
    double term2();
    double correction1();

    double z0();
    double z0R();
    double alpham1();
    double isitz();
    double alpham2(double alpha1);
    double thirdcorrection();
    double secondcorrection();
    double order1fractionalRg2O(const int i, const int j);
    double order2fractionalRg2O(const int i, const int j, const int p, const int q);
    double omatfast(int i, int j, arma::mat &m);  //finds c^T lti c in 1D quickly
    double triomatfast(const int i, const int j, const int u, const int v, arma::mat &m); //finds det(c^T lti c) quickly for 2 delta functions in 1D
    /*     Calculated the summed correction of 3 mayer f functions f_ij f_jk f_ki, which os O(N^4) terms     */
    double threeBondIrreducible(arma::mat &m); 
    arma::mat triomat(const int i, const int j, const int u, const int v, arma::mat &m);
private:
    arma::mat mat2d;
};

// Calculates Laplacian from adjacency matrix
arma::mat calcLaplacian(const arma::mat &adj);
// Calculates the radius of gyration matrix operator
arma::mat rg2mat(int N);
// finds c^T lti c in 1D quickly
arma::mat readMatrix(std::string line, double *dim, double *eps, double *a);
// Create linear matrix adjacency
arma::mat linearAdjacencyMatrix(const int N);
// Create Bond
arma::mat bond(const int i, const int j, const int size);