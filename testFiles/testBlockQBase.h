#include "../src/BlockQBase.h"


void testBlockQBase(){

	DMRGCat::BlockQBase block;
	block.genSiteQBase();
	block.print("single site base");


	DMRGCat::BlockQBase blocks(block);
	DMRGCat::BlockQBase newsys(blocks,block);
	newsys.print("kron base");
	system("pause");
}