#include "SuperBlock.h"
#include "Conjugate.h"


DMRGCat::SuperBlock::SuperBlock(Parameter& para, Block& sys, Block& m, Block& n, Block& env){
	std::vector<int> qno = {para.getL()/2,para.getL()/2};
	TotQNoID = DMRGCat::getID(qno); ;
	PToS = &sys;
	PToM = &m;
	PToN = &n;
	PToE = &env;
	Para = &para;

	Dim = GsWave.setWave(TotQNoID,sys,m,n,env);	
	GsWave0.setWave(TotQNoID, sys, m, n, env);

	calGroundState();
}

void DMRGCat::SuperBlock::f1tof2(const double *f1, double *f2){
	GsWave0.v2QWave(f1);
	in2out(GsWave0, GsWave);
	GsWave.QWave2v(f2);
}

void DMRGCat::SuperBlock::in2out(const QWave& in, QWave& out){
	//SiteH sys,m,n,env
	in.oneBody(BlockS, PToS->QOperator.at(SiteH), out);
	in.oneBody(BlockM, PToM->QOperator.at(SiteH), out);
	in.oneBody(BlockN, PToM->QOperator.at(SiteH), out);
	in.oneBody(BlockE, PToE->QOperator.at(SiteH), out);
	
	//sys-m
	in.twoBody(BlockS, PToS->QOperator.at(CupDag), BlockM, PToM->QOperator.at(Cup), Para->getT(), out);
	in.twoBody(BlockS, PToS->QOperator.at(Cup), BlockM, PToM->QOperator.at(CupDag), Para->getT(), out);
	in.twoBody(BlockS, PToS->QOperator.at(CdownDag), BlockM, PToM->QOperator.at(Cdown), Para->getT(), out);
	in.twoBody(BlockS, PToS->QOperator.at(Cdown), BlockM, PToM->QOperator.at(CdownDag), Para->getT(), out);
	in.twoBody(BlockS, PToS->QOperator.at(Nup), BlockM, PToM->QOperator.at(Ndown), Para->getU(), out);

	//m-n
	in.twoBody(BlockM, PToS->QOperator.at(CupDag), BlockN, PToM->QOperator.at(Cup), Para->getT(), out);
	in.twoBody(BlockM, PToS->QOperator.at(Cup), BlockN, PToM->QOperator.at(CupDag), Para->getT(), out);
	in.twoBody(BlockM, PToS->QOperator.at(CdownDag), BlockN, PToM->QOperator.at(Cdown), Para->getT(), out);
	in.twoBody(BlockM, PToS->QOperator.at(Cdown), BlockN, PToM->QOperator.at(CdownDag), Para->getT(), out);
	in.twoBody(BlockM, PToS->QOperator.at(Nup), BlockN, PToM->QOperator.at(Ndown), Para->getU(), out);

	//n-env
	in.twoBody(BlockN, PToS->QOperator.at(CupDag), BlockE, PToM->QOperator.at(Cup), Para->getT(), out);
	in.twoBody(BlockN, PToS->QOperator.at(Cup), BlockE, PToM->QOperator.at(CupDag), Para->getT(), out);
	in.twoBody(BlockN, PToS->QOperator.at(CdownDag), BlockE, PToM->QOperator.at(Cdown), Para->getT(), out);
	in.twoBody(BlockN, PToS->QOperator.at(Cdown), BlockE, PToM->QOperator.at(CdownDag), Para->getT(), out);
	in.twoBody(BlockN, PToS->QOperator.at(Nup), BlockE, PToM->QOperator.at(Ndown), Para->getU(), out);
}





void DMRGCat::SuperBlock::calGroundState(){
	Conjugate con(Dim);
	con.ErrorBar = 4.e-20;
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
		//con.restart(iter);
		//}
	}

	double GsEnergy = con.eng;
	con.NormTo1(con.f0);
	GsWave.v2QWave(con.f0.memptr());// change sup.Wave = GSWave <- f0
}