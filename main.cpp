// g++ main.cpp multibody.cpp -larmadillo -o multibody
#include <regex>
#include <iostream>
#include <armadillo>
#include <sstream>
#include <vector>
#include "multibody.h"

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
        std::cout << chain.rg20 << " " << res << " " << 0 << " " << std::endl; //chain.thirdcorrection() << std::endl;
        std::cerr << chain.rg20 << " " << res << " " << 0 << " " << std::endl; 
        //std::cerr << "Third Virial " << std::pow(chain.secondcorrection(), 2) - 2.0*chain.thirdcorrection() << std::endl;
        return 0; //Only run once
    }
    return 0;
}