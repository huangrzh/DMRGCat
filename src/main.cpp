#include <ostream>
#include <istream>
#include "../testFiles/testU1Q.h"
#include "../testFiles/testBlockQBase.h"
#include "U1Q.h"
#include "Parameter.h"
#include "Block.h"


int main(){	
	
	std::ifstream ifile("SaveData/Para.txt");
	
	DMRGCat::Parameter para;
	para.load(ifile);
	ifile.close();
	para.print();

	DMRGCat::Block block(para);
	block.print();

#ifdef VISUAL
	system("pause");
#endif
	return 0;
}