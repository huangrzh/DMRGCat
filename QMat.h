#ifndef QMat_H
#define QMat_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <armadillo>
#include <ostream>
#include "Model.h"
#include "BlockQBase.h"

namespace DMRGCat{

class QMat{
	
public:
	QMat();
	QMat(const std::vector<std::pair<int,int>>&, const std::vector<double>&);
		void matCompress(arma::mat& newMat, int dimL, const arma::mat& old, int dimR);
	~QMat();


	QMat(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase);
	friend std::ostream& operator<<(std::ostream&, const QMat&);

	void clear();
	void save(std::ofstream& savefile)const;
	void load(std::ifstream& loadfile);
private:
	std::unordered_map<int, int> RQID2MatNo;
	std::unordered_map<int, int> LQID2MatNo;
	std::unordered_map<int, int> R2LID;
	std::vector<arma::mat> SubMat;
};


}
#endif
