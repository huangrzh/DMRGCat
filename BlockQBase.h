#ifndef BLOCKQBASE_H
#define BLOCKQBASE_H

#include <iostream>
#include <unordered_map>
#include <map>
#include "Model.h"

namespace DMRGCat{

class BlockQBase{
public:
	~BlockQBase();
	BlockQBase();
	
	BlockQBase(const DMRGCat::Model& model);
	BlockQBase(const BlockQBase &b1, const BlockQBase &b2);
	
	
	friend void getKronOrder(const BlockQBase &b1, const BlockQBase &b2, std::map<std::pair<int,int>,int>& startDim);
	void truncate(BlockQBase &ubase);
	friend std::ostream& operator<<(std::ostream& os, const BlockQBase& base);
	friend class QMat;
private:
	std::unordered_map<int, int> SubQIDDim;

	//StartDim is a temp var, and can be cleared for memory reason. 
	//When we need this var again, we can generate it again via getKronOrder any time;
	std::map<std::pair<int,int>, int> StartDim;
};

}
#endif
