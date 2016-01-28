#ifndef MULCHAIN_H
#define MULCHAIN_H
#include <armadillo>
#include <string>
#include "Parameter.h"


namespace DMRGCat{

	class MulChain{

	public:
		MulChain(std::ifstream& infile);
		void constructChain();

		bool IsNeighbor(int, int);


		std::string XEdge;
		std::string YEdge;
		int Lx, Ly;


		int no2X(int);
		int no2Y(int);
		int xy2No(int, int);
		void getNeightbor(int, std::vector<int>&);
		void print();
		void print(std::string s);
	private:
		std::string PathMethod;
		arma::Mat<int> XY2No;
		arma::Col<int> No2X;
		arma::Col<int> No2Y;
	};

}



#endif