#include "SuperBlock.h"



DMRGCat::SuperBlock::SuperBlock(Parameter& para, Block& sys, Block& m, Block& n, Block& env){
	TotQNo = para.getParticleNo();
	PToSys = &sys;
	PToM = &m;
	PToN = &n;
	PToEnv = &env;
	Para = &para;

	Dim = GsWave.setWave(TotQNo,sys,m,n,env);	
}



