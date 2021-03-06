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
	
	void update(const Parameter& para, const Block& old);
	void initial(const Parameter& para);
	void operator=(const Block& var);
	void reNorm(const BlockQBase& UBase, const QMat& reNormU);
	void clear();


	void save(std::ofstream& savefile)const;
	void save(int s)const;
	void load(std::ifstream& loadfile);
	void load(int s);
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
