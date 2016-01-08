#include "../BlockQBase.h"

void testBlockQBase(){
	DMRGCat::BlockQBase block;
	std::cout << "=======================\n";
	std::cout << "test kron\n";
	DMRGCat::BlockQBase blocks(block);
	DMRGCat::BlockQBase newsys(blocks,block);
	std::cout << newsys << std::endl;
}