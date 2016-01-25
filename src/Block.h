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
	Block();
	Block(const Parameter& para);
	Block(const Block& var);
	~Block();
	Block(const Parameter& para, const Block& old);

	void initial(const Parameter& para);
	void operator=(const Block& var);
	void trunc(const BlockQBase& UBase, const QMat& truncU);
	void clear();


	void save(std::ofstream& savefile)const;
	void load(std::ifstream& loadfile);
	void print()const;
	void print(std::string)const;


	friend std::ostream& operator<<(std::ostream&, const Block&);
	friend class QWave;
	friend class SuperBlock;
private:
	BlockQBase QSpace;
    std::vector<QMat> QOperator;
    
};


}
#endif
