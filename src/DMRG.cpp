#include "DMRG.h"
#include "SuperBlock.h"



void DMRGCat::DMRG::initial(){
	std::cout << __FUNCTION__ << ": \n\n";
	std::ifstream ifile("SaveData/Para.txt");
	Para.load(ifile);
	ifile.close();
	Para.print();

	SubS.initial(Para);
	SubM.initial(Para);
	SubN.initial(Para);
	SubE.initial(Para);
}



void DMRGCat::DMRG::warmUp(){
	int lpart = Para.getL() / 2;

	for (int i = 1; i < lpart; i++){
		SuperBlock superChain(Para, SubS, SubM, SubN, SubE);
	}
}