#include "datatype.h"
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


	//CupDag
	std::vector<std::pair<int, int>> idvec;
	idvec.push_back({ id10, id00 });
	idvec.push_back({ id11, id01 });
	QMat cupd(idvec);
	QOperator.push_back(cupd);



	//Cdown
	idvec.clear();
	idvec.push_back({id00,id01});
	idvec.push_back({id10,id11});
	std::vector<double> coe = {1,-1};
	QMat cdown(idvec,coe);
	QOperator.push_back(cdown);


	//CdownDag
	idvec.clear();
	idvec.push_back({ id01, id00 });
	idvec.push_back({ id11, id10 });
	std::vector<double> coe = { 1, -1 };
	QMat cdown(idvec, coe);
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



DMRGCat::Block::Block(const Parameter& para, const Block& old){
	Block added(para);
	QSpace.kron(old.QSpace,added.QSpace);

	QMat temp;
	QOperator = std::vector<QMat>(8,temp);

	//eye
	QOperator.at(Eye).kron(old.QOperator.at(Eye), added.QOperator.at(Eye), QSpace);

	//Up
	QOperator.at(Cup).kron(old.QOperator.at(Eye), added.QOperator.at(Cup), QSpace);

	QOperator.at(CupDag).trans(QOperator.at(Cup));

	//Down
	QOperator.at(Cdown).kron(old.QOperator.at(Eye), added.QOperator.at(Cdown), QSpace);

	QOperator.at(CdownDag).trans(QOperator.at(Cdown));
	
	//Nup
	QOperator.at(Nup).kron(old.QOperator.at(Eye), added.QOperator.at(Nup), QSpace);

	//Ndown
	QOperator.at(Ndown).kron(old.QOperator.at(Eye), added.QOperator.at(Ndown), QSpace);

	//Hamiltonian;
	QOperator.at(SiteH).kron(old.QOperator.at(SiteH), added.QOperator.at(Eye), QSpace);
	DMRGCat::QMat tempO;	
	tempO.kron(old.QOperator.at(Eye), added.QOperator.at(SiteH), QSpace);
	QOperator.at(SiteH).add(tempO,QSpace);

	tempO.kron(old.QOperator.at(CupDag), added.QOperator.at(Cup), QSpace);
	tempO.time(para.getT());
	QOperator.at(SiteH).add(tempO,QSpace);
	tempO.trans();
	QOperator.at(SiteH).add(tempO,QSpace);

	tempO.kron(old.QOperator.at(CdownDag), added.QOperator.at(Cdown), QSpace);
	tempO.time(para.getT());
	QOperator.at(SiteH).add(tempO, QSpace);
	tempO.trans();
	QOperator.at(SiteH).add(tempO, QSpace);

	tempO.kron(old.QOperator.at(Nup), added.QOperator.at(Ndown), QSpace);
	tempO.time(para.getU());
	QOperator.at(SiteH).add(tempO, QSpace);



	????????????????????
}


void DMRGCat::Block::trunc(const BlockQBase& UBase, const QMat& truncU){
	QSpace.truncate(UBase);
	for (auto& x : QOperator){
		x.trunc(UBase, truncU);
	}
}