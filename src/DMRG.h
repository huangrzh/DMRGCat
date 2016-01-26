#ifndef DMRG_H
#define DMRG_H

#include <vector>
#include "Parameter.h"
#include "Block.h"



namespace DMRGCat{

class DMRG{
	public:
		DMRG();
		void initial();
		void warmUp();
		void sweep();		
	private:
		Parameter Para;
		Block SubS;
		Block SubM;
		Block SubN;
		Block SubE;

};

	
}


#endif

