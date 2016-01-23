#include <iostream>
#include "../src/U1Q.h"
#include "../src/QMat.h"


void testQMat(){
	std::cout << "test QMat\n";
	//q(up,down)
	int id00 = DMRGCat::getID({ 0, 0 });
	int id10 = DMRGCat::getID({ 1, 0 });
	int id01 = DMRGCat::getID({ 0, 1 });
	int id11 = DMRGCat::getID({ 1, 1 });
	//Cup
	std::vector<std::pair<int, int>> idvec;
	idvec.push_back({ id00, id10 });
	idvec.push_back({ id01, id11 });
	DMRGCat::QMat cup(idvec);
	cup.print("c up operator");
	system("pause");

	//CupDag
	std::vector<std::pair<int, int>> idvec2;
	idvec2.push_back({ id10, id00 });
	idvec2.push_back({ id11, id01 });
	DMRGCat::QMat cupd(idvec2);
	cupd.print("c up dag");
	system("pause");



	DMRGCat::BlockQBase block;
	block.genSiteQBase();
	DMRGCat::BlockQBase blocks(block);
	DMRGCat::BlockQBase newsys(blocks, block);
	DMRGCat::QMat kronMat;
	kronMat.kron(cupd, cup, newsys);
	kronMat.print("kron mat");
	system("pause");


	std::cout << "\n\n";
}