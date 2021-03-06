#include <iostream>
#include <armadillo>
#include <sstream>
#include <vector>
#include <assert.h>
#include <cmath>
#include "multibody.h"


#define _USE_MATH_DEFINES
 
arma::mat bond(const int i, const int j, const int size){
	arma::mat adj = arma::zeros(size,size);
	adj(i,j) = 1.0;  adj(j,i) = 1.0;
	return calcLaplacian(adj);
}

//Calculates the Laplacian from the Adjacency matrix
arma::mat calcLaplacian(const arma::mat &adj) {
    return arma::diagmat(arma::sum(adj)) - adj;
}
 
//Calculates the radius of gyration matrix operator
arma::mat rg2mat(int N) {
    arma::mat adj = arma::ones(N, N) - arma::eye(N, N);
    return 1.0 / ((double)N*N)*calcLaplacian(adj);
}
 
//finds c^T lti c in 1D quickly
double GaussSystem::omatfast(int i, int j, arma::mat &m) {
    return m.at(i, i) + m.at(j, j) - 2.0*m.at(i, j);
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

/* 
	Create the adjacency matrix for a linear polymer
*/
arma::mat linearAdjacencyMatrix(const int N) {
	arma::vec prototype = arma::zeros(N);
	prototype(1) = 1;
	return arma::toeplitz(prototype);
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
    return rv;
}
 
double GaussSystem::threeBondIrreducible(arma::mat &m) {
	arma::mat rv = arma::zeros(3,3);
	int size = m.n_cols;	//size of matrix, it should be square
	arma::mat c = arma::zeros(size, 3); //c matrix for 3 delta functions
	arma::mat ct = arma::zeros(3, size); //c transpose
	double acc = 0; int dumb;
	for(int i = 0; i < N-1; i++){
		for(int j = i+1; j < N; j++) {
			for(int k = 0; k < N; k++){ if(k != i && k != j) {
				// Stuff
				acc += std::pow(triomatfast(i,j,j,k, m), -D/2.0);
			}}
		}
	}
	return acc;
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

double GaussSystem::order1fractionalRg2O(const int i, const int j) {
	double o = omatfast(i, j, lti);
	double on = omatfast(i, j, opm);
	return 1 - on/o/rg20;
}

double GaussSystem::order2fractionalRg2O(const int i, const int j, const int p, const int q) {
	arma::mat o = triomat(i, j, p, q, lti);
	arma::mat on = triomat(i, j, p, q, opm);
	return 1 - arma::trace(arma::inv(o)*on)/rg20;
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

/*
 Calculate third virial coefficient
 calculates pairs f_{i,j}f_{u,v} where i < j, v > j, u >= i
 a heavy goddamn calculation
*/
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
double GaussSystem::correction1() {
    double acc = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            acc += std::pow(omatfast(i, j, lti), -D/2.0);
        }
    }
    return acc;
}
 
//Calculates \alpha - 1
double GaussSystem::term1() {
    double acc = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            acc += std::pow(omatfast(i, j, lti), -D/2.0)*(1 - order1fractionalRg2O(i,j));
        }
    }
    return acc;
}

double GaussSystem::term2() {
    double acc = 0;
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            for (int u = i; u < N - 1; u++) {
                for (int v = u + 1; v < N; v++) {
                    if (i < u || j < v) {
                        acc += std::pow(triomatfast(i, j, u, v, lti), -D/2.0)*(order2fractionalRg2O(i,j,u,v) - 1);
                    }
                }
            }
        }
    }
    return correction1()*term1() + acc;
}

//Calculates \alpha - 1
double GaussSystem::alpham1() {
    double acc = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            acc += omatfast(i, j, opm)*std::pow(omatfast(i, j, lti), -(1.0*D + 2.0) / 2.0) ;
        }
    }
    return (2*std::pow(2.0*M_PI*a*a / D, -D / 2.0)/rg20)*acc;
}

double GaussSystem::isitz() {
    double acc = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            acc += std::pow(omatfast(i, j, lti), -1.0*D/2.0) ;
        }
    }
    return (2*std::pow(2.0*M_PI*a*a / D, -D / 2.0))*acc;
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
    return term1 - 4*term2/rg20*2.0*M_PI*a*a / D ; 
    
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
 

 
