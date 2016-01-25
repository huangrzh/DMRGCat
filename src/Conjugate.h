/*
Conjugate gradient method for diagonalizing
a sparse symmetric matrix
*/
#include "armadillo"

#ifndef CONJUGATE_H
#define CONJUGATE_H

namespace DMRGCat{

	class Conjugate {
	public:
		Conjugate(const int &Dim);
		~Conjugate();
		double ErrorBar;

		void NormTo1(arma::vec& f); // normalize <f|f> to 1
		double f1timesf2(const arma::vec& f, const arma::vec& g);
		int abc_2(const long &iter);
		void abc_4();

		arma::vec f0;
		arma::vec f1;
		arma::vec f2;
		arma::vec f3;
		long Dim;
		double eng;   // EigenValue of Hamiltonian
		void Ini_f0();
		void restart(long &iter);
		void setErrorBar(double error);
		double getErrorBar(const arma::vec& fin, const arma::vec& fout);
	private:

		double y00;
		double y01;
		double y02;
		double y22;
		double x33_old;
	};
}

#endif
