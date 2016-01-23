#include "U1Q.h"


DMRGCat::U1Q::U1Q(){
#ifdef TWO_Q
  Num = 2;
  ID = MAX_Q*Q.at(0) + Q.at(1);
#endif
#ifdef ONE_Q
  Num = 1;
  ID = Q.at(0);
#endif
  Q.resize(Num);
}


DMRGCat::U1Q::U1Q(int inputid){
	ID = inputid;
#ifdef TWO_Q
	Num = 2;
	Q.resize(Num);
	int max = MAX_Q;
	Q[0] = ID/max;
	Q[1] = ID%max;
#endif
#ifdef ONE_Q
	Num = 1;
	Q.resize(Num);
	Q[0] = ID;
#endif
	
}


DMRGCat::U1Q::U1Q(const std::vector<int>& inputQ){
#ifdef TWO_Q
	Num = 2;
	Q.resize(Num);
	Q = inputQ;
	ID = MAX_Q*Q.at(0) + Q.at(1);
#endif
#ifdef ONE_Q
	Num = 1;
	Q.resize(Num);
	Q = inputQ;
	ID = Q.at(0);
#endif
	try {
		if (!(inputQ.size() == Num)) {
			throw std::runtime_error("Error in U1Q(std::vector<int>: inputQ.size!=Num)");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
	
}


DMRGCat::U1Q::U1Q(int* inputQ){	
#ifdef TWO_Q
	Num = 2;
	Q.resize(2);
	Q[0] = inputQ[0];
	Q[1] = inputQ[1];
	ID = MAX_Q*inputQ[0] + inputQ[1];
#endif
#ifdef ONE_Q
	Num = 1;
	Q.resize(1);
	Q[0] = inputQ[0];
	ID = inputQ[0];
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


int DMRGCat::U1Q::getChargeNo()const{
	int charge = 0;
	for (const auto& x : Q){
		charge += x;
	}
	return charge;
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




int DMRGCat::getID(const std::vector<int>& var){
#ifdef TWO_Q
	return MAX_Q*var.at(0) + var.at(1);
#endif
#ifdef ONE_Q
	return var.at(0);
#endif
}



#ifdef FERMION


bool DMRGCat::hasSign(int id){
#ifdef TWO_Q
	int q1 = id / MAX_Q;
	int q2 = id%MAX_Q;
	if ((q1 + q2) % 2 == 0){
		return false;
	}
	else{
		return true;
	}
#endif
#ifdef ONE_Q
	if (id % 2 == 1){
		return true;
	}
	else{
		return false;
	}
#endif
}



bool DMRGCat::hasSign(int lqid, int rqid){
#ifdef TWO_Q
	int lq1 = lqid / MAX_Q;
	int lq2 = lqid % MAX_Q;
	int rq1 = rqid / MAX_Q;
	int rq2 = rqid % MAX_Q;
	int deltaq = lq1 + lq2 - rq1 - rq2;
	if (deltaq < 0){
		deltaq = -deltaq;
	}
	if (deltaq % 2 == 1){
		return true;
	}
	else{
		return false;
	}
#endif
#ifdef ONE_Q
	if ((lqid + rqid) % 2 == 1){
		return true;
	}
	else{
		return false;
	}
#endif
}
int DMRGCat::getFermionSign(int id){
#ifdef TWO_Q
	int q1 = id / MAX_Q;
	int q2 = id%MAX_Q;
	if ((q1 + q2) % 2 == 0){
		return 1;
	}
	else{
		return -1;
	}
#endif
#ifdef ONE_Q
	if (id % 2 == 1){
		return -1;
	}
	else{
		return 1;
	}
#endif
}

int DMRGCat::getFermionSign(int lqid, int rqid){
#ifdef TWO_Q
	int lq1 = lqid / MAX_Q;
	int lq2 = lqid % MAX_Q;
	int rq1 = rqid / MAX_Q;
	int rq2 = rqid % MAX_Q;
	if ((lq1 + rq1) % 2 == 1 || (lq2 + rq2) % 2 == 1){
		return -1;
	}
	else{
		return 1;
	}
#endif
#ifdef ONE_Q
	if ((lqid+rqid) % 2 == 1){
		return -1;
	}
	else{
		return 1;
	}
#endif
}


#endif



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
#ifdef ONE_Q
       output<<std::setw(3)<<GQvar.Q[0]<<")   ";
#endif
#ifdef TWO_Q
       for(int i=0;i<1;i++){
          output << GQvar.Q[i] <<" , ";
       }
       output << GQvar.Q[1] <<")  ";
#endif
    return output;
}