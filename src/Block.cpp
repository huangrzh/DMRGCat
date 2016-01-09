#include "Block.h"
#include "setting.h"


//Single site block;
DMRGCat::Block::Block(const Parameter& para){
#ifdef FERMION_HUBBARD
	//q(up,down)
	int id00 = DMRGCat::getID({ 0, 0 });
	int id10 = DMRGCat::getID({ 1, 0 });
	int id01 = DMRGCat::getID({ 0, 1 });
	int id11 = DMRGCat::getID({ 1, 1 });
	//eye
	QMat eyeOp;
	eyeOp.eyeQMat(QSpace);
	QOperator.push_back(eyeOp);

	//Cup
	std::vector<std::pair<int, int>> idvec;
	idvec.push_back({id00,id10});
	idvec.push_back({id01,id11});
	QMat cup(idvec);
	QOperator.push_back(cup);

	//Cdown
	idvec.clear();
	idvec.push_back({id00,id01});
	idvec.push_back({id10,id11});
	std::vector<double> coe = {1,-1};
	QMat cdown(idvec,coe);
	QOperator.push_back(cdown);


	//Nup
	idvec.clear();
	idvec.push_back({id10,id10});
	idvec.push_back({id11,id11});
	QMat Nup(idvec);
	QOperator.push_back(Nup);


	//Ndown
	idvec.clear();
	idvec.push_back({ id01, id01 });
	idvec.push_back({ id11, id11 });
	QMat Ndown(idvec);
	QOperator.push_back(Ndown);


	//Site Hamiltonian
	idvec.clear();
	idvec.push_back({id11,id11});
	coe.clear();
	coe.push_back(para.getU());
	QMat siteH(idvec,coe);
	QOperator.push_back(siteH);
#endif
}



DMRGCat::Block::Block(const Block& old, const Block& added){
	QSpace.kron(old.QSpace,added.QSpace);

	QMat temp;
	QOperator = std::vector<QMat>(6,temp);

	//eye
	QOperator.at(0).kron(old.QOperator.at(0), added.QOperator.at(0), QSpace);

	//Up, Down, N, Hamiltonian;

}


void DMRGCat::Block::trunc(const BlockQBase& UBase, const QMat& truncU){
	QSpace.truncate(UBase);
	for (auto& x : QOperator){
		x.trunc(UBase, truncU);
	}
}