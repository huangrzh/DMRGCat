#include <iostream>
#include <fstream>
#include <functional>
#include "../src/QWave.h"


void testQWave(){
	std::cout << __FUNCTION__ << ": \n\n";

	std::ifstream ifile("SaveData/Para.txt");
	DMRGCat::Parameter para;
	para.load(ifile);
	ifile.close();
	para.print();

	DMRGCat::Block sys(para);
	//sys.print("sys block");
	//system("pause");

	DMRGCat::Block m(sys);
	DMRGCat::Block n(sys);
	DMRGCat::Block env(sys);

	std::vector<int> qno = { para.getL() / 2, para.getL() / 2 };
	int qnoid = DMRGCat::getID(qno);

	DMRGCat::QWave psi(qnoid, sys, m, n, env);
	std::cout << "qwave constructor over\n";
	system("pause");
	psi.print("psi");
}