#include "DMRG.h"
#include "SuperBlock.h"


DMRGCat::DMRG::DMRG(){
	std::cout << __FUNCTION__ << ": \n\n";
	initial();
	warmUp();
	sweep();
}
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
	std::cout << __FUNCTION__ << ": \n\n";
	int lpart = Para.getL() / 2;

	int Os = 1;
	int Oe = Para.getL();
	for (int i = 1; i < lpart; i++){
		SuperBlock superChain(Para, SubS, SubM, SubN, SubE);
		
		QMat reNormU;
		BlockQBase UBase;
		superChain.getReNormU(reNormU, UBase);
		
		Block olds(SubS);
		SubS.update(Para,SubS);
		SubS.reNorm(UBase, reNormU);
		SubS.save(Os);
		SubE.save(Oe);

		Os++;
		Oe--;
	}
}



void DMRGCat::DMRG::sweep(){
	std::cout << __FUNCTION__ << ": \n\n";
	int SweepDir = 1;
	int Os = Para.getL() / 2 - 1;
	for (int i = 1; i <= Para.getSweepNo(); i++){
		bool OneSweep = { true };		
		while (OneSweep){
			SuperBlock superChain(Para, SubS, SubM, SubN, SubE);

			QMat reNormU;
			BlockQBase UBase;
			superChain.getReNormU(reNormU, UBase);

			Block olds(SubS);
			SubS.update(Para, SubS);
			SubS.reNorm(UBase, reNormU);
			SubS.save(Os);

			Os++;

			OneSweep = (SweepDir == 1 && Os == Para.getL() - 3) || (SweepDir == -1 && Os == 3);
		}
		Os = SweepDir == 1 ? Para.getL() : 1;
		SweepDir = -SweepDir;
	}
}