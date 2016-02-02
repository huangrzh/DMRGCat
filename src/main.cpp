#include <ostream>
#include <istream>
#include "MulChain.h"
#include "DMRG.h"
#include "Block2D.h"


int main(){	
	std::ifstream loadMulChain("SaveData/MulChainData.txt", std::ios::in);
	DMRGCat::MulChain chain0(loadMulChain);
	loadMulChain.close();

	chain0.print("chain0");
	system("pause");
	DMRGCat::DMRG task0;

#ifdef VISUAL
	system("pause");
#endif
	return 0;
}