#ifndef UONEQNO_H
#define UONEQNO_H


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "setting.h"

namespace DMRGCat{


//for fermion hubbard modle, (up,down)
class U1Q{
	
public:
//construct U1Q
    U1Q();
	U1Q(int);
	U1Q(const std::vector<int>& inputQ);
	U1Q(int*);
    ~U1Q();
	
	
//U1Q access
	int getNoOfQ()const;
	int getQ(int)const;
	int getID()const;
	int getChargeNo()const;
	
	
//operatoration of U1Q
    U1Q& operator= (const U1Q&);
    U1Q& operator+=(const U1Q&);
    U1Q& operator-=(const U1Q&);
    U1Q  operator+ (const U1Q&)const;
    U1Q  operator- (const U1Q&)const;
    bool operator==(const U1Q&)const;
    bool operator!=(const U1Q&)const;
    bool operator> (const U1Q&)const;
    bool operator< (const U1Q&)const; 
   
	friend std::ostream& operator<<(std::ostream&, const U1Q&);

private:
    std::vector<int> Q;
    int Num;
	int ID;
};

#ifdef FERMION
bool hasSign(int qid);
bool hasSign(int lid, int rid);
int getFermionSign(int lqid, int rqid);
int getFermionSign(int id);
#endif

int getID(const std::vector<int>&);
int getAddID(int id1, int id2);

}
#endif
