#include <time.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <omp.h>
using namespace std;
#include "Conjugate.h"
DMRGCat::Conjugate::~Conjugate(){}


DMRGCat::Conjugate::Conjugate(const int &dim) {
	Dim = dim;
	f0.zeros(Dim);
	f1.zeros(Dim);
	f2.zeros(Dim);
	f3.zeros(Dim);
	Ini_f0();

};



void DMRGCat::Conjugate::restart(long &iter){
	iter = 0;
	NormTo1(f0);
}

/*
void Conjugate::Start(char y_n) {
CreateSpace() ;
if ( y_n == 'y' ) Ini_f0() ; // randomly initialize f0 if y_n=='y'
}
*/

void DMRGCat::Conjugate::Ini_f0() {
	srand(time(NULL));
#pragma omp parallel for
	for (long i = 0; i < Dim; i++){
		f0[i] = (double)rand();
	}
}



double DMRGCat::Conjugate::getErrorBar(const arma::vec& fin, const arma::vec& fout){
	double inner1 = f1timesf2(fin, fin);   // y00 = <f0|f0>
	double inner2 = f1timesf2(fin, fout);	// y01 = <f0|H|f0>
	double lamda = inner2 / inner1;    // energy	

	arma::vec newv(Dim);
	double aux = 1.0 / inner1;
#pragma omp parallel for
	for (long j = 0; j < Dim; j++){
		newv[j] = (fout[j] - lamda * fin[j]) / inner1;
	}

	double x33 = f1timesf2(newv, newv);
	return (inner1*x33) / (lamda*lamda);
}



int DMRGCat::Conjugate::abc_2(const long &iter) {
	y00 = f1timesf2(f0, f0);   // y00 = <f0|f0>
	y01 = f1timesf2(f0, f1);	// y01 = <f0|H|f0>
	eng = y01 / y00;    // energy	

	double aux = 2.0 / y00;
#pragma omp parallel for
	for (long j = 0; j < Dim; j++){
		f3[j] = (f1[j] - eng * f0[j]) * aux;
	}
	double x33 = f1timesf2(f3, f3); // x33=<f3|f3>
	double errorNow = (y00 * x33) / (eng * eng);
	if (iter % 30 == 0){
		cout << "    Energy[" << iter << "] = " << setprecision(18) << eng << ", errorbar = " << errorNow << endl;
	}

	if (errorNow < ErrorBar) return 1;
	if (iter == 0) {
#pragma omp parallel for
		for (long i = 0; i < Dim; i++){
			f2[i] = -f3[i];
		}
	}
	else {
		double aaa = x33 / x33_old;
#pragma omp parallel for
		for (long k = 0; k < Dim; k++){
			f2[k] = -f3[k] + aaa * f2[k];
		}
	}
	x33_old = x33;
	y02 = f1timesf2(f0, f2); // y02=<f0|f2>
	y22 = f1timesf2(f2, f2);  // y22=<f2|f2>
	return 0;
}

void DMRGCat::Conjugate::abc_4() {
	double x03 = f1timesf2(f0, f3);
	double x23 = f1timesf2(f2, f3);
	double xa = x23 * y02 - x03 * y22;
	double xb = x23 * y00 - y01 * y22;
	double xc = x03 * y00 - y01 * y02;
	double alpha = (-xb + sqrt(xb * xb - 4.0 * xa * xc)) / (2.0 * xa);

#pragma omp parallel for
	for (long i = 0; i < Dim; i++) {
		f0[i] += alpha * f2[i];
		f1[i] += alpha * f3[i];
		f3[i] = 0.0;
	}
}

double DMRGCat::Conjugate::f1timesf2(const arma::vec& f, const arma::vec& g) {
	double x = 0.0;
#pragma omp parallel for reduction(+:x)
	for (long i = 0; i < Dim; i++){
		x += f[i] * g[i];
	}
	return x;
}




void DMRGCat::Conjugate::NormTo1(arma::vec& f) {
	double x = sqrt(f1timesf2(f, f));
#pragma omp parallel for
	for (long i = 0; i < Dim; i++){
		f[i] /= x;
	}
}



void DMRGCat::Conjugate::setErrorBar(double error){
	if (error < 4.e-28){
		ErrorBar = 4.e-28;
	}
	else{
		ErrorBar = error;
	}
}
