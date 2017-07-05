#include <armadillo>

// Class that handles perturbation theory of a polymer
class GaussSystem {
public:
    int N;
    arma::mat lap, lti, rg2, opm;
    double rg20, D, a;
 
    double _z0 = -1.0;
 
    GaussSystem(arma::mat laplacian, double dimension, double seglen);
    double z0();
    double z0R();
    double alpham1();
    double alpham2(double alpha1);
    double thirdcorrection();
    double secondcorrection();
 
    double omatfast(int i, int j, arma::mat &m);  //finds c^T lti c in 1D quickly
    double triomatfast(const int i, const int j, const int u, const int v, arma::mat &m); //finds det(c^T lti c) quickly for 2 delta functions in 1D
    arma::mat triomat(const int i, const int j, const int u, const int v, arma::mat &m);
private:
    arma::mat mat2d;
};

arma::mat calcLaplacian(const arma::mat &adj);
// Calculates the radius of gyration matrix operator
arma::mat rg2mat(int N);
// finds c^T lti c in 1D quickly
arma::mat readMatrix(std::string line, double *dim, double *eps, double *a);
// Create linear matrix adjacency
arma::mat linearAdjacencyMatrix(const int N);
// Create Bond
arma::mat bond(const int i, const int j, const int size);