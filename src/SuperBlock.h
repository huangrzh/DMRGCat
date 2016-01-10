#ifndef SUPERBLOCK_H
#define SUPERBLOCK_H

#include <vector>
#include "Parameter.h"
#include "QMat.h"
#include "Block.h"
#include "QWave.h"



namespace DMRGCat{

class SuperBlock{
	public:
		SuperBlock(Parameter& para, Block& sys, Block& m, Block& n, Block& env);
		
	private:
		QWave GsWave;
		int TotQNo;
		int Dim;
		Block* PToSys;
		Block* PToM;
		Block* PToN;
		Block* PToEnv;
		Parameter* Para;		
		
		
		void calGroundState();
		void calPhysicalQuantity();
		
		void preWaveStep1();
		void preWaveStep2();
};

	
}


#endif

