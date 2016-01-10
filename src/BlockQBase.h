#ifndef BLOCKQBASE_H
#define BLOCKQBASE_H

#include <iostream>
//#include <unordered_map>
#include <map>

namespace DMRGCat{

class BlockQBase{
public:
	~BlockQBase();
	BlockQBase();
	void genSiteQBase();
	BlockQBase(const BlockQBase &b1, const BlockQBase &b2);
	void kron(const BlockQBase &b1, const BlockQBase &b2);



	void save(std::ofstream& savefile)const;
	void load(std::ifstream& loadfile);
	
	void truncate(const BlockQBase &ubase);
	friend std::ostream& operator<<(std::ostream& os, const BlockQBase& base);
	friend void getKronOrder(const BlockQBase &b1, const BlockQBase &b2, std::map<std::pair<int, int>, int>& startDim);
	friend class QMat;
	friend class Block;
	friend class QWave;
private:
	//We set SubQIDDim to make it in defined order. Thus we won't have chaos when update site;
	std::map<int, int> SubQIDDim;

	//StartDim is a temp var, and can be cleared for memory reason. 
	//When we need this var again, we can generate it again via getKronOrder any time;
	std::map<std::pair<int,int>, int> StartDim;
};

}
#endif
