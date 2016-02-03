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

		std::string xEdge()const{
			return XEdge;
		}

		std::string yEdge()const{
			return YEdge;
		}

		
		int getLx()const;
		int getLy()const;

		int no2X(int)const;
		int no2Y(int)const;
		int xy2No(int, int)const;
		void getNeightbor(int, std::vector<int>&)const;
		void print();
		void print(std::string s);
	private:
		std::string XEdge;
		std::string YEdge;
		int Lx, Ly;
		std::string PathMethod;
		arma::Mat<int> XY2No;
		arma::Col<int> No2X;
		arma::Col<int> No2Y;
	};


	

}



#endif