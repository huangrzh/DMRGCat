#include "SuperBlock.h"
#include "Conjugate.h"


DMRGCat::SuperBlock::SuperBlock(Parameter& para, Block& sys, Block& m, Block& n, Block& env){
	TotQNoID = para.getTotQNoID();
	PToS = &sys;
	PToM = &m;
	PToN = &n;
	PToE = &env;
	Para = &para;

	Dim = GsWave.setWave(TotQNoID,sys,m,n,env);	
	GsWave0.setWave(TotQNoID, sys, m, n, env);
	calGroundState();
}



void DMRGCat::SuperBlock::calGroundState(){
	Conjugate con(Dim);
	con.setErrorBar(4.e-20);
	long iter = 0;
	bool breakfor = false;
	for (int j = 0; j < 400; j++){
		if (iter == 0){
			f1tof2(con.f0.memptr(), con.f1.memptr());// f1 = H f0
		}
		breakfor = con.abc_2(iter);
		if (breakfor){
			break;
		}	
		f1tof2(con.f2.memptr(), con.f3.memptr()); // f3 = H f2
		con.abc_4();
		iter++;

		//if (iter == 30){
		//	con.restart(iter);
		//}
	}

	double GsEnergy = con.eng;
	con.NormTo1(con.f0);
	GsWave.v2QWave(con.f0.memptr());// change sup.Wave = GSWave <- f0
	std::cout << ", Energy = " << con.eng << "\n";
}



void DMRGCat::SuperBlock::f1tof2(const double *f1, double *f2){	
	GsWave.zeros();
	GsWave0.v2QWave(f1);
	in2out(GsWave0, GsWave);
	GsWave.QWave2v(f2);
}



void DMRGCat::SuperBlock::in2out(const QWave& in, QWave& out){
	//SiteH sys,m,n,env
	in.oneBody(BlockS, PToS->QOperator.at(SiteH), out);	
	in.oneBody(BlockM, PToM->QOperator.at(SiteH), out);
	in.oneBody(BlockN, PToN->QOperator.at(SiteH), out);
	in.oneBody(BlockE, PToE->QOperator.at(SiteH), out);
	
	
	//sys-m
	in.twoBody(BlockS, PToS->QOperator.at(CupDag),   BlockM, PToM->QOperator.at(Cup),       Para->getT(),  out);
	in.twoBody(BlockS, PToS->QOperator.at(Cup),		 BlockM, PToM->QOperator.at(CupDag),   -Para->getT(),  out);
	in.twoBody(BlockS, PToS->QOperator.at(CdownDag), BlockM, PToM->QOperator.at(Cdown),     Para->getT(),  out);
	in.twoBody(BlockS, PToS->QOperator.at(Cdown),	 BlockM, PToM->QOperator.at(CdownDag), -Para->getT(),  out);
	

	
	//m-n
	in.twoBody(BlockM, PToM->QOperator.at(CupDag),   BlockN, PToN->QOperator.at(Cup),		Para->getT(), out);
	in.twoBody(BlockM, PToM->QOperator.at(Cup),      BlockN, PToN->QOperator.at(CupDag),   -Para->getT(), out);
	in.twoBody(BlockM, PToM->QOperator.at(CdownDag), BlockN, PToN->QOperator.at(Cdown),     Para->getT(), out);
	in.twoBody(BlockM, PToM->QOperator.at(Cdown),    BlockN, PToN->QOperator.at(CdownDag), -Para->getT(), out);
	
	
	
	//n-env
	in.twoBody(BlockN, PToN->QOperator.at(CupDag),   BlockE, PToE->QOperator.at(Cup),       Para->getT(), out);
	in.twoBody(BlockN, PToN->QOperator.at(Cup),      BlockE, PToE->QOperator.at(CupDag),   -Para->getT(), out);
	in.twoBody(BlockN, PToN->QOperator.at(CdownDag), BlockE, PToE->QOperator.at(Cdown),     Para->getT(), out);
	in.twoBody(BlockN, PToN->QOperator.at(Cdown),    BlockE, PToE->QOperator.at(CdownDag), -Para->getT(), out);
}


void DMRGCat::SuperBlock::getReNormU(QMat& U, BlockQBase& UBase)const{
	QMat waveMat;
	GsWave.wave2QMat(PToS->QSpace, PToM->QSpace, PToN->QSpace, PToE->QSpace, waveMat);
	waveMat.getReNormUAndBase(Para->getD(), U, UBase);
}