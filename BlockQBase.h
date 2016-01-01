#ifndef BLOCKQBASE_H
#define BLOCKQBASE_H

#include "U1Q.h"

namespace DMRGCat{

class BlockQBase{
	BlockQBase();
	~BlockQBase();
	
	void kron(const BlockQBase &b1,const BlockQBase &b2);
	void truncate(BlockQBase &ubase);
	friend std::ostream& operator<<(std::ostream& os, const BlockQBase& base);
private:
	int NoOfQ;
	std::vector<U1Q> BlockQ;
	std::vector<int> QDim;
};

}
#endif
