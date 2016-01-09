#ifndef SUPERBLOCK_H
#define SUPERBLOCK_H

#include <vector>
#include "Parameter.h"
#include "QMat.h"
#include "Block.h"



namespace DMRGCat{

class SuperBlock{
	public:
		SuperBlock(const Parameter& para, const Block& sys, const Block& m, const Block& n, const Block& env);
		
		void getTruncU(BlockQBase& UBase, QMat& truncU);
	private:
		Block* PToSys;
		Block* PToM;
		Block* PToN;
		Block* PToEnv;

		std::vector<std::pair<int,int>> MNQID;
		std::vector<QMat> SysEnvQMat;
		
		
		void calGroundState();
		void calPhysicalQuantity();
		
		void preWaveStep1();
		void preWaveStep2();
};

	
}


#endif

