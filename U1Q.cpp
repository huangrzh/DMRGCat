#include "U1Q.h"

DMRGCat::U1Q::U1Q(){
  Num=NUMBER_OF_GOOD_QUANTUM_NUMBER;
  Q.resize(Num);
#ifdef TWO_Q
  ID = MAX_Q*Q.at(0) + Q.at(1);
#endif
#ifdef ONE_Q
  ID = Q.at(0);
#endif
}

DMRGCat::U1Q::U1Q(int inputid){
	Num = NUMBER_OF_GOOD_QUANTUM_NUMBER;
	Q.resize(Num);
	ID = inputid;
#ifdef TWO_Q
	int max = MAX_Q;
	Q[0] = ID/max;
	Q[1] = ID%max;
#endif
#ifdef ONE_Q
	Q[0] = ID;
#endif
}


DMRGCat::U1Q::U1Q(const std::vector<int>& inputQ){
	Num = NUMBER_OF_GOOD_QUANTUM_NUMBER;
	try {
		if (!(inputQ.size()==Num)) { 
			throw std::runtime_error("Error in U1Q(std::vector<int>: inputQ.size!=Num)");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}	
	Q.resize(Num);
	Q = inputQ;
#ifdef TWO_Q
	ID = MAX_Q*Q.at(0) + Q.at(1);
#endif
#ifdef ONE_Q
	ID = Q.at(0);
#endif
}


DMRGCat::U1Q::U1Q(int* inputQ){
	Num = NUMBER_OF_GOOD_QUANTUM_NUMBER;
	std::vector<int> qs(inputQ, inputQ + Num);
	Q = qs;
#ifdef TWO_Q
	ID = MAX_Q*Q.at(0) + Q.at(1);
#endif
#ifdef ONE_Q
	ID = Q.at(0);
#endif
}



DMRGCat::U1Q::~U1Q(){}


int DMRGCat::U1Q::getNoOfQ()const{
	return Num;
}


int DMRGCat::U1Q::getQ(int i)const{
	return Q.at(i);
}

int DMRGCat::U1Q::getID()const{
	return ID;
}


int DMRGCat::getAddID(int id1, int id2){
#ifdef TWO_Q
	int q1 = id1/MAX_Q + id2/MAX_Q;
	int q2 = id1%MAX_Q + id2%MAX_Q;
	return MAX_Q*q1 + q2;
#endif
#ifdef ONE_Q
	return id1 + id2;
#endif
}

DMRGCat::U1Q& DMRGCat::U1Q::operator= (const DMRGCat::U1Q& Gvar){
	Q = Gvar.Q;
	ID = Gvar.ID;
    return *this;
}


bool DMRGCat::U1Q::operator==(const DMRGCat::U1Q& Gvar)const{
    if(ID==Gvar.ID) return false;
    return true;
}


bool DMRGCat::U1Q::operator!=(const DMRGCat::U1Q& Gvar)const{
	if (ID != Gvar.ID) return false;
	return true;
}


bool DMRGCat::U1Q::operator> (const DMRGCat::U1Q& Gvar)const{
    if(ID>Gvar.ID) return true ;     
    return false;
}


bool DMRGCat::U1Q::operator< (const DMRGCat::U1Q& Gvar)const{
	if (ID<Gvar.ID) return true;
	return false;
}



DMRGCat::U1Q& DMRGCat::U1Q::operator+=(const DMRGCat::U1Q& Gvar){
    for(int i=0;i<Num;i++){
       Q[i]+=Gvar.Q[i];
    }
#ifdef TWO_Q
	ID = MAX_Q*Q.at(0) + Q.at(1);
#endif
#ifdef ONE_Q
	ID = Q.at(0);
#endif
    return *this;
}



DMRGCat::U1Q& DMRGCat::U1Q::operator-=(const DMRGCat::U1Q& Gvar){
  for(int i=0;i<Num;i++){
    Q[i]-=Gvar.Q[i];
  }
#ifdef TWO_Q
  ID = MAX_Q*Q.at(0) + Q.at(1);
#endif
#ifdef ONE_Q
  ID = Q.at(0);
#endif
  return *this;
}

DMRGCat::U1Q  DMRGCat::U1Q::operator+ (const DMRGCat::U1Q& Gvar)const{
    DMRGCat::U1Q Gsum=*this;
    Gsum+=Gvar;
    return Gsum;
}


DMRGCat::U1Q  DMRGCat::U1Q::operator- (const DMRGCat::U1Q& Gvar)const{
  DMRGCat::U1Q Gminus=*this;
  Gminus-=Gvar;
  return Gminus;
}

std::ostream& DMRGCat::operator<<(std::ostream& output,const DMRGCat::U1Q& GQvar){
    output<< "(";
    if(NUMBER_OF_GOOD_QUANTUM_NUMBER==1){
       output<<std::setw(3)<<GQvar.Q[0]<<")   ";
    } else{
       for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER-1;i++){
          output<<std::setw(3)<<GQvar.Q[i]<<",";
       }
       output<<std::setw(3)<<GQvar.Q[NUMBER_OF_GOOD_QUANTUM_NUMBER-1]<<"  )";
    }
    return output;
}