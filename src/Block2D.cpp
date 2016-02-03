#include "datatype.h"
#include "Block2D.h"


DMRGCat::Block2D::Block2D(const Block2D& var) : Block(var){
	SiteNo = var.SiteNo;
	NoOfSites = var.NoOfSites;
	InSite2OperatorNo = var.InSite2OperatorNo;
	OutSite2OperatorNo = var.OutSite2OperatorNo;
	OutQOCd = var.OutQOCd;
	OutQOCdDag = var.OutQOCdDag;
	OutQOCu = var.OutQOCu;
	OutQOCuDag = var.OutQOCuDag;
	InQOCd = var.InQOCd;
	InQOCdDag = var.InQOCdDag;
	InQOCu = var.InQOCu;
	InQOCuDag = var.InQOCuDag;
}


void DMRGCat::Block2D::operator=(const Block2D& var){
	clear2D();
	QSpace = var.QSpace;
	QOperator = var.QOperator;
	SiteNo = var.SiteNo;
	NoOfSites = var.NoOfSites;
	InSite2OperatorNo = var.InSite2OperatorNo;
	OutSite2OperatorNo = var.OutSite2OperatorNo;
	OutQOCd = var.OutQOCd;
	OutQOCdDag = var.OutQOCdDag;
	OutQOCu = var.OutQOCu;
	OutQOCuDag = var.OutQOCuDag;
	InQOCd = var.InQOCd;
	InQOCdDag = var.InQOCdDag;
	InQOCu = var.InQOCu;
	InQOCuDag = var.InQOCuDag;
}


void DMRGCat::Block2D::initial2D(const Parameter& para, int site){
	initial(para);
	SiteNo = site;
	NoOfSites = 1;
}


DMRGCat::Block2D::Block2D(const Parameter& para, int site){
	initial2D(para, site);
}


void DMRGCat::Block2D::update2D(const Parameter& para, const MulChain& chain, const Block2D& added, const Block2D& old){
	clear2D();
	SiteNo = added.SiteNo;
	NoOfSites = old.NoOfSites + 1;
	update(para, old);
	if (old.OutSite2OperatorNo.size() == chain.getLy()){
		int sweep_dir = added.SiteNo > old.SiteNo ? 1 : -1;
		int sitex_ = chain.no2X(added.SiteNo);
		int sitex = sitex_ - sweep_dir;
		int sitey = chain.no2Y(added.SiteNo);
		int leftSiteNo = chain.xy2No(sitex, sitey);
		int left_OP_no = old.OutSite2OperatorNo.at(leftSiteNo);
		QMat tempO;
		tempO.kron(added.QOperator.at(CupDag), old.QOperator.at(Cup), QSpace);
		tempO.time(para.getT());
		QOperator.at(SiteH).add(tempO, QSpace);
		tempO.trans();
		QOperator.at(SiteH).add(tempO, QSpace);
	}

		
	OutSite2OperatorNo = old.OutSite2OperatorNo;
	OutQOCd = old.OutQOCd;
	OutQOCdDag = old.OutQOCdDag;
	OutQOCu = old.OutQOCu;
	OutQOCuDag = old.OutQOCuDag;

	InSite2OperatorNo = old.InSite2OperatorNo;
	InQOCd = old.InQOCd;
	InQOCdDag = old.InQOCdDag;
	InQOCu = old.InQOCu;
	InQOCuDag = old.InQOCuDag;


	//Out size == 0
	if (old.OutSite2OperatorNo.size() == 0){
		OutSite2OperatorNo[old.SiteNo] = 0;		
		OutQOCd.push_back(added.QOperator.at(Cdown));
		OutQOCdDag.push_back(added.QOperator.at(CdownDag));
		OutQOCu.push_back(added.QOperator.at(Cup));
		OutQOCuDag.push_back(added.QOperator.at(CupDag));
	}


	int oldsize = OutSite2OperatorNo.size();
	QMat tempO0;
	for (int i = 0; i < oldsize; i++){
		tempO0.kron(added.QOperator.at(Eye), old.OutQOCd.at(i), QSpace);
		OutQOCd.at(i) = tempO0;
		OutQOCdDag.at(i).trans(OutQOCd.at(i));
		tempO0.kron(added.QOperator.at(Eye), old.OutQOCu.at(Cup), QSpace);
		OutQOCu.at(i) = tempO0;
		OutQOCuDag.at(i).trans(OutQOCu.at(i));
	}

	//Out size >0, <Ly
	if (old.OutSite2OperatorNo.size() < chain.getLy()){		
		OutSite2OperatorNo[SiteNo] = oldsize + 1;
		OutQOCd.push_back(added.QOperator.at(Cdown));
		OutQOCdDag.push_back(added.QOperator.at(CdownDag));
		OutQOCu.push_back(added.QOperator.at(Cup));
		OutQOCuDag.push_back(added.QOperator.at(CupDag));
	}

	//Out size == Ly
	else{		
		int old_in_size = old.InSite2OperatorNo.size();
		QMat tempO1;
		for (int i = 0; i < old_in_size; i++){
			tempO1.kron(added.QOperator.at(Eye), old.InQOCd.at(i), QSpace);
			InQOCd.at(i) = tempO1;
			InQOCdDag.at(i).trans(InQOCd.at(i));
			tempO1.kron(added.QOperator.at(Cup), old.QOperator.at(i), QSpace);
			InQOCu.at(i) = tempO1;
			InQOCuDag.at(i).trans(InQOCd.at(i));
		}
		//old.In size == 0 || old.In size == Ly
		if (old_in_size == 0 || old_in_size<chain.getLy() ){
			InSite2OperatorNo[SiteNo] = old_in_size;			
			InQOCd.push_back(added.QOperator.at(Cdown));
			InQOCdDag.push_back(added.QOperator.at(CdownDag));
			InQOCu.push_back(added.QOperator.at(Cup));
			InQOCuDag.push_back(added.QOperator.at(CupDag));
		}
		else{
			int sweepDir = SiteNo>old.SiteNo ? 1 : -1;
			int sitex_ = chain.no2X(SiteNo);
			int sitey = chain.no2Y(SiteNo);
			int sitex = sitex_ - sweepDir;
			int old_new_no = chain.xy2No(sitex, sitey);
			int OPNo = InSite2OperatorNo(old_new_no);
			InQOCd.at(OPNo) = added.QOperator.at(Cdown);
			InQOCdDag.at(OPNo) = added.QOperator.at(CdownDag);
			InQOCu.at(OPNo) = added.QOperator.at(Cup);
			InQOCuDag.at(OPNo) = added.QOperator.at(CupDag);
			InSite2OperatorNo.erase(old_new_no);
			InSite2OperatorNo[SiteNo] = OPNo;
		}
	}
}



void DMRGCat::Block2D::clear2D(){
	clear();
	InSite2OperatorNo.clear();
	OutSite2OperatorNo.clear();
	OutQOCd.clear();
	OutQOCdDag.clear();
	OutQOCu.clear();
	OutQOCuDag.clear();
	InQOCd.clear();
	InQOCdDag.clear();
	InQOCu.clear();
	InQOCuDag.clear();
}



void DMRGCat::Block2D::reNorm2D(const BlockQBase& UBase, const QMat& reNormU){
	reNorm(UBase, reNormU);

	int osize = OutSite2OperatorNo.size();
	int isize = InSite2OperatorNo.size();
	for (int i = 0; i < osize; i++){
		OutQOCd.at(i).trunc(UBase, reNormU);
		OutQOCdDag.at(i).trunc(UBase, reNormU);
		OutQOCu.at(i).trunc(UBase, reNormU);
		OutQOCuDag.at(i).trunc(UBase, reNormU);
	}
	for (int i = 0; i < isize; i++){
		InQOCd.at(i).trunc(UBase, reNormU);
		InQOCdDag.at(i).trunc(UBase, reNormU);
		InQOCu.at(i).trunc(UBase, reNormU);
		InQOCuDag.at(i).trunc(UBase, reNormU);
	}
}


void DMRGCat::Block2D::save2D(std::ofstream& savefile)const{
	save(savefile);
	
	savefile.write((const char*)&NoOfSites, sizeof(int));
	int osize = OutSite2OperatorNo.size();
	int isize = InSite2OperatorNo.size();
	savefile.write((const char*)&osize, sizeof(int));
	for (const auto& x : OutSite2OperatorNo){
		savefile.write((const char*)&(x.first), sizeof(int));
		savefile.write((const char*)&(x.second), sizeof(int));
	}

	savefile.write((const char*)&isize, sizeof(int));
	for (const auto& x : InSite2OperatorNo){
		savefile.write((const char*)&(x.first), sizeof(int));
		savefile.write((const char*)&(x.second), sizeof(int));
	}

	for (int i = 0; i < osize; i++){
		OutQOCd.at(i).save(savefile);
		OutQOCu.at(i).save(savefile);			
	}
	for (int i = 0; i < isize; i++){
		InQOCd.at(i).save(savefile);
		InQOCu.at(i).save(savefile);
	}
}

void DMRGCat::Block2D::save2D()const{
	std::string outs = "Data/";
	outs += std::to_string(SiteNo);
	std::ofstream sfout(outs, std::ios::binary | std::ios::out);
	save2D(sfout);
	sfout.close();
}



void DMRGCat::Block2D::load2D(std::ifstream& loadfile){
	clear2D();
	load(loadfile);

	loadfile.read((char*)&NoOfSites, sizeof(int));
	unsigned int osize = 0;
	unsigned int isize = 0;	
	loadfile.read((char*)&osize, sizeof(int));
	for (int i = 0; i < osize; i++){
		int xfirst = 0;
		int xsecond = 0;
		loadfile.read((char*)&xfirst, sizeof(int));
		loadfile.read((char*)&xsecond, sizeof(int));
		OutSite2OperatorNo[xfirst] = xsecond;
	}

	loadfile.read((char*)&isize, sizeof(int));
	for (int i = 0; i < isize; i++){
		int xfirst = 0;
		int xsecond = 0;
		loadfile.read((char*)&xfirst, sizeof(int));
		loadfile.read((char*)&xsecond, sizeof(int));
		InSite2OperatorNo[xfirst] = xsecond;
	}

	QMat tempO;
	OutQOCd = std::vector<QMat>(osize, tempO);
	OutQOCdDag = std::vector<QMat>(osize, tempO);
	OutQOCu = std::vector<QMat>(osize, tempO);
	OutQOCuDag = std::vector<QMat>(osize, tempO);
	for (int i = 0; i < osize; i++){
		OutQOCd.at(i).load(loadfile);
		OutQOCu.at(i).load(loadfile);
		OutQOCdDag.at(i).trans(OutQOCd.at(i));
		OutQOCuDag.at(i).trans(OutQOCu.at(i));
	}

	InQOCd = std::vector<QMat>(isize, tempO);
	InQOCdDag = std::vector<QMat>(isize, tempO);
	InQOCu = std::vector<QMat>(isize, tempO);
	InQOCuDag = std::vector<QMat>(isize, tempO);
	for (int i = 0; i < isize; i++){
		InQOCd.at(i).load(loadfile);
		InQOCu.at(i).load(loadfile);
		InQOCdDag.at(i).trans(InQOCd.at(i));
		InQOCuDag.at(i).trans(InQOCu.at(i));
	}
}


void DMRGCat::Block2D::load2D(int s){
	std::string ins = "Data/";
	ins += std::to_string(s);
	std::ifstream sfin(ins, std::ios::binary | std::ios::in);
	load2D(sfin);
	sfin.close();
}


void DMRGCat::Block2D::print2D()const{
	std::cout << "site number = " << SiteNo << "\n";
	std::cout << "No of sites = " << NoOfSites << "\n";
	print();
	std::cout << "Out size = " << OutSite2OperatorNo.size() << "\n";
	std::cout << "In size = " << InSite2OperatorNo.size() << "\n";
}


void DMRGCat::Block2D::print2D(std::string s)const{
	std::cout << s << "\n";
	print2D();
}