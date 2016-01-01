#include <string>
#include <vector>
#include <iostream>
#include "../Model.h"
#include "../U1Q.h"


void testModel(){
	std::string modelname = "FermionHubbard";
	DMRGCat::Model FHubbard(modelname);
	std::vector<int> qids;
	FHubbard.getSpaceQID(qids);
	for (const auto& x : qids){
		DMRGCat::U1Q q(x);
		std::cout << q << std::endl;
	}
}