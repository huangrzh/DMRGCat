#ifndef BLOCK_H
#define BLOCK_H


#include <vector>
#include "Parameter.h"
#include "U1Q.h"
#include "BlockQBase.h"
#include "QMat.h"

namespace DMRGCat{


class Block{
	
public:
	Block(const Parameter& para);
	~Block();
	Block(const Block& old, const Block& added);


	void trunc(const BlockQBase& UBase, const QMat& truncU);
	void clear();


	void save(std::ofstream& savefile)const;
	void load(std::ifstream& loadfile);
	friend std::ostream& operator<<(std::ostream&, const Block&);
	friend class QWave;
private:
	BlockQBase QSpace;
    std::vector<QMat> QOperator;
    
};


}
#endif
