#include "Parameter.h"


DMRGCat::Parameter::Parameter(){}



void DMRGCat::Parameter::load(std::ifstream& s){
#ifdef FERMION_HUBBARD
	char x = 'q';
	while (x != '=') s >> x;
	s >> T;


	x = 'q';
	while (x != '=') s >> x;
	s >> U;


	x = 'q';
	while (x != '=') s >> x;
	s >> ParticleNo;



	x = 'q';
	while (x != '=') s >> x;
	s >> SavedD;

	x = 'q';
	while (x != '=') s >> x;
	s >> SweepNo;
#endif
}



std::ostream& DMRGCat::operator<<(std::ostream& s, const DMRGCat::Parameter& para){
#ifdef FERMION_HUBBARD
    s << "Parameter: \n\n" 


      << "t = " << para.T << std::endl 

      << "u = " << para.U << std::endl

      << "SavedD = " << para.SavedD  << std::endl 
	  
	  << "ParticleNo = " << para.ParticleNo  << std::endl 
 
      << "MaxSweepNo = " << para.SweepNo << std::endl;
#endif
	      
    return s;
}



double DMRGCat::Parameter::getU()const{
	return U;
}
double DMRGCat::Parameter::getT()const{
	return T;
}
int DMRGCat::Parameter::getD()const{
	return SavedD;
}
int DMRGCat::Parameter::getSweepNo()const{
	return SweepNo;
}
int DMRGCat::Parameter::getParticleNo()const{
	return ParticleNo;
}









