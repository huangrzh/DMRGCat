#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>
#include <ostream>
#include <fstream>
#include "setting.h"
//class Parameter is dependant on model studied.

namespace DMRGCat{

class Parameter{
	public:
		Parameter();
		void print();
		void load(std::ifstream& loadfile);
		friend std::ostream& operator<<(std::ostream& s, const Parameter& para);

		void setTotQNoID(int id);
		int getL()const;
		int getD()const;
		int getTotQNoID()const;
		int getSweepNo()const;
#ifdef FERMION_HUBBARD
		double getU()const;
		double getT()const;		
#endif
	private:

#ifdef FERMION_HUBBARD
		int ChainL;
		double T;
		double U;
		int TotQNoID;
		int SavedD;
		int SweepNo; 		// no need to store;
#endif
		
};

	
}


#endif

