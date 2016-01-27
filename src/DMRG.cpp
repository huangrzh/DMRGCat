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
	std::cout << "\n\n";
}



void DMRGCat::DMRG::warmUp(){
	std::cout << __FUNCTION__ << ": \n\n";
	int lpart = Para.getL() / 2;

	int Os = 1;
	int Oe = Para.getL();
	SubS.save(Os);
	SubE.save(Oe);

	for (int i = 1; i < lpart-1; i++){
		std::cout << "@" << std::setiosflags(std::ios::left) << std::setw(4) << Os;
		std::cout << "@" << std::setiosflags(std::ios::left) << std::setw(4) << Oe;
		std::vector<int> qno = {Os+1,Os+1};
		Para.setTotQNoID(DMRGCat::getID(qno));
		SuperBlock superChain(Para, SubS, SubM, SubN, SubE);
		
		QMat reNormU;
		BlockQBase UBase;
		superChain.getReNormU(reNormU, UBase);
		
		Block olds(SubS);
		SubS.update(Para,olds);
		SubS.reNorm(UBase, reNormU);
		SubE = SubS;
		
		Os++;
		Oe--;
		SubS.save(Os);
		SubE.save(Oe);
	}
	std::cout << "\n\n";
}



void DMRGCat::DMRG::sweep(){
	std::cout << __FUNCTION__ << ": \n\n";
	int SweepDir = 1;
	int Os = Para.getL() / 2 - 1;
	std::vector<int> qno = { Os + 1, Os + 1 };
	Para.setTotQNoID(DMRGCat::getID(qno));
	for (int i = 1; i <= Para.getSweepNo(); i++){
		bool OneSweep = { true };		
		while (OneSweep){
			std::cout << "Sweep#" << std::setiosflags(std::ios::left) << std::setw(4) << i;
			std::cout << "@" << std::setiosflags(std::ios::left) << std::setw(4) << Os;
			std::cout << "@" << std::setiosflags(std::ios::left) << std::setw(4) << Os + 3 * SweepDir;

			SuperBlock superChain(Para, SubS, SubM, SubN, SubE);

			QMat reNormU;
			BlockQBase UBase;
			superChain.getReNormU(reNormU, UBase);
			Block olds(SubS);
			SubS.update(Para, olds);
			SubS.reNorm(UBase, reNormU);
			
			Os += SweepDir;
			SubS.save(Os);
			SubE.load(Os + 3 * SweepDir);
			
			OneSweep = !((SweepDir == 1 && Os == Para.getL() - 3) || (SweepDir == -1 && Os == 4));
		}
		std::cout << "\n\n";
		Os = SweepDir == 1 ? Para.getL() : 1;
		SweepDir = -SweepDir;
		Block tempBlock;
		std::swap(SubS, SubE);
	}
}