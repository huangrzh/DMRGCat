#ifndef QMat_H
#define QMat_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <armadillo>
#include <ostream>
#include "setting.h"
#include "BlockQBase.h"

namespace DMRGCat{

class QMat{
	
public:
	QMat();                        //LRIDs
	QMat(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<double>&);
	QMat(const std::vector<std::pair<int, int>>& LRIDs);
	QMat(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<std::pair<int,int>>& LRDims);
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

#ifdef FERMION
	bool getIsFermion()const;
#endif

	friend class QWave;


	//-----------------------------------------------------------------------------------------
	// --------------------------operations----------------------------------------------------
	// out += lamda in
	friend void time(const double& lamda, const QMat& in, QMat& out);

	// out += O * in
	friend void leftTime(const QMat& O, const QMat& in, QMat& out);

	// out += lamda O * in
	friend void leftTime(const double& lamda, const QMat& O, const QMat& in, QMat& out);

	// out += in * O
	friend void rightTime(const QMat& O, const QMat& in, QMat& out);

	// out += lamda in * O
	friend void rightTime(const double& lamda, const QMat& O, const QMat& in, QMat& out);


	// out += leftO * in * rightO
	friend void lrTime(const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out);

	// out += lamda leftO * in * rightO
	friend void lrTime(const double& lamda, const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out);
	// --------------------------operations----------------------------------------------------
	//-----------------------------------------------------------------------------------------
	

private:

#ifdef FERMION
	bool IsFermion;
#endif

	std::unordered_map<int, int> RQID2MatNo;
	std::unordered_map<int, int> LQID2MatNo;
	std::unordered_map<int, int> R2LID;
	std::vector<std::pair<int, int>> LRID;
	std::vector<arma::mat> SubMat;
};


}
#endif
