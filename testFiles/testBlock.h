#include <iostream>
#include <fstream>
#include <functional>
#include "../src/Block.h"


void testBlock(){
	std::cout << __FUNCTION__ << ": \n\n";

	std::ifstream ifile("SaveData/Para.txt");
	DMRGCat::Parameter para;
	para.load(ifile);
	ifile.close();
	para.print();

	DMRGCat::Block sys(para);
	sys.print("sys");
	system("pause");

	DMRGCat::Block newSys(para, sys);
	newSys.print("newSys");

	std::ofstream ofile("Data/testSaveBlock",std::ios::out);
	newSys.save(ofile);
	ofile.close();

	system("pause");

	DMRGCat::Block iblock;
	std::ifstream ifile_block("Data/testSaveBlock");
	iblock.load(ifile_block);
	ifile.close();
	iblock.print("iblock");
}