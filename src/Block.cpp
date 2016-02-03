#include "datatype.h"
#include "Block.h"
#include "setting.h"


//Single site block;
DMRGCat::Block::Block(const Parameter& para){
	initial(para);
}


void DMRGCat::Block::initial(const Parameter& para){
	clear();

	QSpace.genSiteQBase();
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
	idvec.push_back({ id00, id10 });
	idvec.push_back({ id01, id11 });
	QMat cup(idvec);
	QOperator.push_back(cup);


	//CupDag
	QMat cupd;
	QOperator.push_back(cupd);
	QOperator.back().trans(cup);



	//Cdown
	idvec.clear();
	idvec.push_back({ id00, id01 });
	idvec.push_back({ id10, id11 });
	std::vector<double> coe = { 1, -1 };
	QMat cdown(idvec, coe);
	QOperator.push_back(cdown);


	//CdownDag
	QMat cdown2;
	QOperator.push_back(cdown2);
	QOperator.back().trans(cdown);


	//Nup
	idvec.clear();
	idvec.push_back({ id10, id10 });
	idvec.push_back({ id11, id11 });
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
	idvec.push_back({ id11, id11 });
	coe.clear();
	coe.push_back(para.getU());
	QMat siteH(idvec, coe);
	QOperator.push_back(siteH);
#endif
}


DMRGCat::Block::Block(const Block& var){
	QSpace = var.QSpace;
	QOperator = var.QOperator;
}


void DMRGCat::Block::operator=(const Block& var){
	clear();
	QSpace = var.QSpace;
	QOperator = var.QOperator;
}

DMRGCat::Block::Block(const Parameter& para, const Block& old){
	update(para, old);
}


void DMRGCat::Block::update(const Parameter& para, const Block& old){
	clear();
	
	Block added(para);
	QSpace.kron(added.QSpace, old.QSpace);

	QMat temp;
	QOperator = std::vector<QMat>(8, temp);

	//eye
	QOperator.at(Eye).kron(added.QOperator.at(Eye), old.QOperator.at(Eye), QSpace);

	//Up
	QOperator.at(Cup).kron(added.QOperator.at(Cup), old.QOperator.at(Eye), QSpace);
	QOperator.at(CupDag).trans(QOperator.at(Cup));

	//Down
	QOperator.at(Cdown).kron(added.QOperator.at(Cdown), old.QOperator.at(Eye), QSpace);
	QOperator.at(CdownDag).trans(QOperator.at(Cdown));

	//Nup
	QOperator.at(Nup).kron(added.QOperator.at(Nup), old.QOperator.at(Eye), QSpace);

	//Ndown
	QOperator.at(Ndown).kron(added.QOperator.at(Ndown), old.QOperator.at(Eye), QSpace);


	//Hamiltonian;
	QOperator.at(SiteH).kron(added.QOperator.at(Eye), old.QOperator.at(SiteH), QSpace);

	DMRGCat::QMat tempO;
	tempO.kron(added.QOperator.at(SiteH), old.QOperator.at(Eye), QSpace);
	QOperator.at(SiteH).add(tempO, QSpace);


	tempO.kron(added.QOperator.at(CupDag), old.QOperator.at(Cup), QSpace);
	tempO.time(para.getT());
	QOperator.at(SiteH).add(tempO, QSpace);

	tempO.trans();
	QOperator.at(SiteH).add(tempO, QSpace);

	tempO.kron(added.QOperator.at(CdownDag), old.QOperator.at(Cdown), QSpace);
	tempO.time(para.getT());
	QOperator.at(SiteH).add(tempO, QSpace);
	tempO.trans();
	QOperator.at(SiteH).add(tempO, QSpace);
}


void DMRGCat::Block::reNorm(const BlockQBase& UBase, const QMat& reNormU){
	QSpace.truncate(UBase);
	for (auto& x : QOperator){
		x.trunc(UBase, reNormU);
	}
}



void DMRGCat::Block::clear(){
	QSpace.clear();
	QOperator.clear();
}

void DMRGCat::Block::print()const{	
	
	for (int i = 0; i < QOperator.size(); i++){
		std::cout << "QONo = " << i << std::endl;
		QOperator.at(i).print();
		system("pause");
	}
}


void DMRGCat::Block::print(std::string s)const{
	std::cout << s << "\n";
	print();
}



void DMRGCat::Block::save(std::ofstream& savefile)const{
	QSpace.save(savefile);
	unsigned int NumOfBlock = QOperator.size();
	savefile.write((const char*)&NumOfBlock, sizeof(int));

	for (const auto& x : QOperator){
		x.save(savefile);
	}
}

void DMRGCat::Block::save(int s)const{
	std::string outs = "Data/";
	outs += std::to_string(s);
	std::ofstream sfout(outs, std::ios::binary | std::ios::out);
	save(sfout);
	sfout.close();
}

void DMRGCat::Block::load(std::ifstream& loadfile){
	clear();
	QSpace.load(loadfile);
	unsigned int size = 0;
	loadfile.read((char*)&size, sizeof(int));
	QMat tempVar;
	QOperator = std::vector<QMat>(size, tempVar);
	for (int i = 0; i < size; i++){
		QOperator.at(i).load(loadfile);
	}
}


void DMRGCat::Block::load(int s){
	std::string ins = "Data/";
	ins += std::to_string(s);
	std::ifstream sfin(ins, std::ios::binary | std::ios::in);
	load(sfin);
	sfin.close();
}