#include <iostream>
#include <armadillo>
#include <sstream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <regex>

#define _USE_MATH_DEFINES
 
//Calculates the Laplacian from the Adjacency matrix
arma::mat calcLaplacian(const arma::mat &adj) {
    return arma::diagmat(arma::sum(adj)) - adj;
}
 
//Calculates the radius of gyration matrix operator
arma::mat rg2mat(int N) {
    arma::mat adj = arma::ones(N, N) - arma::eye(N, N);
    return 1.0 / ((double)N*N)*calcLaplacian(adj);
}
 
//Class that handles perturbation theory of a polymer
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
 
//finds c^T lti c in 1D quickly
double GaussSystem::omatfast(int i, int j, arma::mat &m) {
    return m.at(i, i) + m.at(j, j) - 2.0*m.at(i, j);
}
 
double GaussSystem::triomatfast(const int i, const int j, const int u, const int v, arma::mat &m) {
    //matrix assosciated with \delta(r_i - r_j) \delta(r_u - r_v)
    //matrix is a reference to a pre-allocated matrix.  This function will be called a butt-ton of times
    //so it should be fast if possible
    mat2d(0, 0) = m(i, i) + m(j, j) - m(i, j) - m(j, i); mat2d(0, 1) = m(i, u) + m(j, v) - m(i, v) - m(j, u);
    mat2d(1, 0) = m(u, i) + m(v, j) - m(v, i) - m(u, j); mat2d(1, 1) = m(u, u) + m(v, v) - m(u, v) - m(v, u);
   
    /*
    manual verification
    arma::mat cm = arma::zeros(10, 2);
    cm(i, 0) = 1; cm(j, 0) = -1;
    cm(u, 1) = 1; cm(v, 1) = -1;
    std::cout << cm.t()*m*cm;
    std::cout << mat2d << std::endl; */
    return arma::det(mat2d);
}

arma::mat GaussSystem::triomat(const int i, const int j, const int u, const int v, arma::mat &m) {
    //matrix assosciated with \delta(r_i - r_j) \delta(r_u - r_v)
    //matrix is a reference to a pre-allocated matrix.  This function will be called a butt-ton of times
    //so it should be fast if possible
    arma::mat rv(2,2);
    rv(0, 0) = m(i, i) + m(j, j) - m(i, j) - m(j, i); rv(0, 1) = m(i, u) + m(j, v) - m(i, v) - m(j, u);
    rv(1, 0) = m(u, i) + m(v, j) - m(v, i) - m(u, j); rv(1, 1) = m(u, u) + m(v, v) - m(u, v) - m(v, u);
   
    /*
    manual verification
    arma::mat cm = arma::zeros(10, 2);
    cm(i, 0) = 1; cm(j, 0) = -1;
    cm(u, 1) = 1; cm(v, 1) = -1;
    std::cout << cm.t()*m*cm;
    std::cout << mat2d << std::endl; */
    return rv;
}
 
 
GaussSystem::GaussSystem(arma::mat laplacian, double dimension, double seglen) {
    N = laplacian.n_cols;   //Matrix size.  Assume it is square
    D = dimension;
    a = seglen;
    lap = laplacian;   
    lti = arma::inv(lap);
    rg2 = rg2mat(N);
    opm = lti*rg2*lti;
    rg20 = arma::trace(lti*rg2);
    mat2d = arma::zeros(2, 2);  //2x2 matrix for use in triomatfast
}
 
double GaussSystem::secondcorrection() {
	double acc = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            acc += std::pow(omatfast(i, j, lti), -D/2.0);
        }
    }
    return acc * std::pow(D/(2*M_PI*a*a), D/2.0);
}

//Calculate third virial coefficient
//calculates pairs f_{i,j}f_{u,v} where i < j, v > j, u >= i
//a heavy goddamn calculation
double GaussSystem::thirdcorrection() {
    double sum = 0;
    double factor = std::pow(D / (2*M_PI*a*a), D);
    double det;
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            for (int u = i; u < N - 1; u++) {
                for (int v = u + 1; v < N; v++) {
                    if (i < u || j < v) {
                        det = triomatfast(i, j, u, v, lti);
                        sum += std::pow(det, -D/2.0);
                    }
                }
            }
        }
    }
    return factor*sum;
}
 
//Calculates \alpha - 1
double GaussSystem::alpham1() {
    double acc = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            acc += std::pow(omatfast(i, j, lti), -(D + 2) / 2.0) * omatfast(i, j, opm);
        }
    }
    return (2*std::pow(2.0*M_PI*a*a / D, -D / 2.0)/rg20) * acc;
}

//next order in \alpha
//requires previous order correction as parameter
double GaussSystem::alpham2(double alpha1) {
	double term1 = 4*alpha1*secondcorrection();
	double term2 = 0;
    arma::mat ok(2,2);	//matrices needed for 3rd order calculation
    arma::mat oki(2,2);
    arma::mat okp(2,2);
    //Horrid calculation
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            for (int u = i; u < N - 1; u++) {
                for (int v = u + 1; v < N; v++) {
                    if (i < u || j < v) {
                        ok = triomat(i, j, u, v, lti); //c^T L^{-1} c
                        oki = arma::inv(ok);	//invert it
                        okp = triomat(i, j, u, v, opm); //c^T L^{-1} R_G^2 L^{-1}
                        term2 += std::pow(arma::det(ok), -D/2.0)*arma::trace(oki*okp);
                    }
                }
            }
        }
    }
    return term1 - 0*4*term2/rg20*2.0*M_PI*a*a / D;
}
 
/*
Calculate <R_G^2>_0
*/
double GaussSystem::z0R() {
    if (_z0 < 0) {  //calculate z0 if have not done so already
        _z0 = z0();
    }
    return _z0*rg20;
}
 
/*
z0 returns the non-interacting partition function
*/
double GaussSystem::z0() {
    double num = std::pow(2.0*M_PI*a*a / D, (double)N);
    double den = arma::det(lap);
    _z0 = std::pow(num / den, D / 2.0);
    return _z0;
}
 
arma::mat readMatrix(std::string line, double *dim, double *eps, double *a) {
        std::stringstream ss;
        ss.str(line);
        int width;
        ss >> width; ss >> *dim; ss >> *eps >> *a;
        arma::mat adjacency(width, width);
        int i = 0;
        while(1) {  //read in numbers
            double elem;
            ss >> elem;
            if (ss) {
                adjacency[i] = elem;
                i++;
            }
            else {
                break;
            }
        }
        assert(i == width*width);
        return adjacency;
}
 
int main(int argc, char** argv)
{
    double eps = 1.0;
    double dim = 3.0;
    double a = 1.0;
    std::regex re("^\\w*\\n*$"); //find an empty line with a new line at end
    /*
    read in file format: N dim eps a {an N x N adjacency matrix as list of numbers, no brackets}
    */
    for (std::string line; std::getline(std::cin, line); ) { //process lines
        if(std::regex_match(line, re)) { //on blank line, exit
            return 0;
        }
        arma::mat adjacency = readMatrix(line, &dim, &eps, &a);
        int N = adjacency.n_cols;
        arma::mat delta = eps*arma::eye(N, N); //This anchors all monomers
        if(eps == 0.0){	//Anchor only 1 monomer.  This has the effect of no confining potential
        	delta = 0*delta;
        	delta(0, 0) = 1.0; //This only anchors one monomer
        }
        arma::mat lap = calcLaplacian(adjacency) + delta;
        GaussSystem chain = GaussSystem(lap, dim, a);
        double res = chain.alpham1();
        std::cout << chain.rg20 << " " << res << " " << chain.alpham2(res) << " " << std::endl; //chain.thirdcorrection() << std::endl;
        //std::cerr << "Third Virial " << std::pow(chain.secondcorrection(), 2) - 2.0*chain.thirdcorrection() << std::endl;
        return 0; //Only run once
    }
    return 0;
}