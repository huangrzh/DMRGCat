#ifndef QMat_H
#define QMat_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <armadillo>
#include <ostream>
#include "BlockQBase.h"

namespace DMRGCat{

class QMat{
	
public:
	QMat();                        //LRIDs
	QMat(const std::vector<std::pair<int,int>>&, const std::vector<double>&);
	QMat(const std::vector<std::pair<int, int>>&);
	~QMat();

	void kron(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase);
	QMat(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase);
		void matCompress(arma::mat& newMat, int dimL, const arma::mat& old, int dimR);
	friend std::ostream& operator<<(std::ostream&, const QMat&);


	void eyeQMat(const BlockQBase&);
	void clear();
	void save(std::ofstream& savefile)const;
	void load(std::ifstream& loadfile);


	void trunc(const BlockQBase& UBase, const QMat& truncU);
private:
	std::unordered_map<int, int> RQID2MatNo;
	std::unordered_map<int, int> LQID2MatNo;
	std::unordered_map<int, int> R2LID;
	std::vector<std::pair<int, int>> LRID;
	std::vector<arma::mat> SubMat;
};


}
#endif
