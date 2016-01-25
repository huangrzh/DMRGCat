#include <iostream>
#include <fstream>
#include <functional>
#include "../src/SuperBlock.h"


void testSuperBlock(){
	std::cout << __FUNCTION__ << ": \n\n";

	std::ifstream ifile("SaveData/Para.txt");
	DMRGCat::Parameter para;
	para.load(ifile);
	ifile.close();
	para.print();

	
	DMRGCat::Block sys(para);
	DMRGCat::Block m(sys);
	DMRGCat::Block n(sys);
	DMRGCat::Block env(sys);


	DMRGCat::Block sys2(para, sys);
	DMRGCat::Block sys3(para, sys2);
	DMRGCat::Block env3(sys3);
	DMRGCat::SuperBlock superBlock8(para, sys3, m, n, env3);	
	std::cout << "SuperBlock over\n";
}