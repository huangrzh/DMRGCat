#ifndef DMRG
#define DMRG

#include <vector>
#include "Parameter.h"
#include "QMat.h"
#include "Block.h"



namespace DMRGCat{

class DMRG{
	public:
		void initial();
		void warmUp();
		void sweep();		
	private:
		

};

	
}


#endif

