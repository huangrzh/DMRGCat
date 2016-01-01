#include "../Model.h"
#include "../BlockQBase.h"

void testBlockQBase(){
	std::string s = "FermionHubbard";
	DMRGCat::Model model(s);

	std::cout << "test single site block base\n";
	DMRGCat::BlockQBase block(model);
	std::cout << block << std::endl;

	std::cout << "=======================\n";
	std::cout << "test kron\n";
	DMRGCat::BlockQBase blocks(block);
	DMRGCat::BlockQBase newsys(blocks,block);
	std::cout << newsys << std::endl;
}