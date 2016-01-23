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
		QWave GsWave, GsWave0;
		int TotQNoID;
		int Dim;
		Block *PToS;
		Block *PToM;
		Block *PToN;
		Block *PToE;
		Parameter *Para;		
		
		
		void calGroundState();
		void in2out(const QWave& in, QWave& out);
		void f1tof2(const double *f1, double *f2);
		void calPhysicalQuantity();
		
		void preWaveStep1();
		void preWaveStep2();

};

	
}


#endif

