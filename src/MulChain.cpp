#include "MulChain.h"
#include <math.h>
#include <iomanip>


//*
//*  ^ ->  ^
//*  |  |  |
//*     >--|
DMRGCat::MulChain::MulChain(std::ifstream& s){
	char x = 'q';
	while (x != '=') s >> x;
	s >> XEdge;

	x = 'q';
	while (x != '=') s >> x;
	s >> YEdge;

	x = 'q';
	while (x != '=') s >> x;
	s >> Lx;


	x = 'q';
	while (x != '=') s >> x;
	s >> Ly;

	x = 'q';
	while (x != '=') s >> x;
	s >> PathMethod;

	constructChain();
}


void DMRGCat::MulChain::constructChain(){	
	XY2No.zeros(Lx+2,Ly+2);
	No2X.zeros(Lx*Ly);
	No2Y.zeros(Lx*Ly);

	std::string s_std = "std";
	if (PathMethod == s_std){
		for (int i = 1; i <= Lx; i++){
			for (int j = 1; j <= Ly; j++){
				int no = i % 2 == 1 ? ((i - 1)*Ly + j) : ((i - 1)*Ly + Ly - j + 1);
				XY2No(i, j) = no;
				No2X(no - 1) = i;
				No2Y(no - 1) = j;
			}
		}
	}

	for(int j=0; j<Ly+2; j++){
		XY2No(0,j) = -1;
		XY2No(Lx+1,j) = -1;
	}

	for(int i=0; i<Lx+2; i++){
		XY2No(i,0) = -1;
		XY2No(i,Ly+1) = -1;
	}
}


int DMRGCat::MulChain::no2X(int no){
	return No2X(no-1);
}

int DMRGCat::MulChain::no2Y(int no){
	return No2Y(no-1);
}


int DMRGCat::MulChain::xy2No(int x, int y){
	return XY2No(x,y);
}

bool DMRGCat::MulChain::IsNeighbor(int no1, int no2){
	bool tempbool = no1>=1 && no1<=Lx*Ly && no2>=1 && no2<=Lx*Ly;
	try{
		if (!tempbool){
			throw std::runtime_error("Error in MulChain no1 or no2 out of range");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}


	int x1 = no2X(no1);
	int y1 = no2Y(no1);
	int x2 = no2X(no2);
	int y2 = no2Y(no2);

	if((x1-x2)*(y1-y2)==0){
		if(abs(x1-x2)==1 || abs(y1-y2)==1){
			return true;
		}
	
		if((XEdge=="P" && abs(x1-x2)==Lx-1 && x1!=x2) || (YEdge=="P" && abs(y1-y2)==Ly-1 && y1!=y2)){
			return true;
		}
	}
	return false;
}


void DMRGCat::MulChain::getNeightbor(int no,std::vector<int>& neighbor){
	neighbor.clear();

	if(no<1 || no>Lx*Ly){
		return;
	}
	int x = no2X(no);
	int y = no2Y(no);
	int yup = y+1;
	int ydown = y-1;
	int xup = x+1;
	int xdown = x-1;

	if(yup>Ly && YEdge=="P" && y!=1 && ydown!=1){
		yup = 1;
	}
	if(ydown<1 && YEdge=="P" && y!=Ly && yup!=Ly){
		ydown = Ly;
	}
	if(xup>Lx && XEdge=="P" && x!=1 &&xdown!=1){
		xup = 1;
	}
	if(xdown<1 && XEdge=="P" && x!=Lx && xup!=Lx){
		xdown = Lx;
	}
	
	int no2;
	no2 = xy2No(x,yup);  
	if(no2>0) neighbor.push_back(no2);
	no2 = xy2No(x,ydown); 
	if (no2>0) neighbor.push_back(no2);
	no2 = xy2No(xup,y);   
	if (no2>0) neighbor.push_back(no2);
	no2 = xy2No(xdown,y); 
	if (no2>0) neighbor.push_back(no2);
}



void DMRGCat::MulChain::print(){
	std::cout << "XEdge = " << XEdge << ", YEdge = " << YEdge << "\n";
	std::cout << "Lx = " << Lx << ", Ly = " << Ly << "\n";
	std::cout << "PathMethod = " << PathMethod << "\n";	

	unsigned int size = Lx*Ly;
	for (int i = 1; i <= size; i++){
		std::cout << "@" << std::setiosflags(std::ios::left) << std::setw(5) << i;
		std::cout << std::setiosflags(std::ios::left) << std::setw(3) << no2X(i);
		std::cout << std::setiosflags(std::ios::left) << std::setw(3) << no2Y(i) << "\n";
	}
}

void DMRGCat::MulChain::print(std::string s){
	std::cout << __FUNCTION__ << "\n";
	std::cout << s << "\n";
	print();
}