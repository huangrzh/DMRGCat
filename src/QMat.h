#ifndef QMat_H
#define QMat_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <armadillo>
#include <ostream>
#include <string>
#include "setting.h"
#include "BlockQBase.h"

namespace DMRGCat{

class QMat{
	
public:
	QMat();                       
	QMat(const QMat& var);
	QMat(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<double>&);
	QMat(const std::vector<std::pair<int, int>>& LRIDs);
	QMat(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<std::pair<int,int>>& LRDims);
		void zero(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<std::pair<int, int>>& LRDims);
	~QMat();

	void kron(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase);
	QMat(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase);
		void matCompress(arma::mat& newMat, int dimL, const arma::mat& old, int dimR);
	friend std::ostream& operator<<(std::ostream&, const QMat&);

	void operator=(const QMat& var);
	void eyeQMat(const BlockQBase&);
	void clear();
	void save(std::ofstream& savefile)const;
	void load(std::ifstream& loadfile);


	void trunc(const BlockQBase& UBase, const QMat& truncU);
	void trans(const QMat&);
	void trans();
	void add(const QMat& added, const BlockQBase& space);
	void time(double lamda);
#ifdef FERMION
	bool getIsFermion()const;
#endif

	friend class QWave;
	friend class Block;


	

	//-----------------------------------------------------------------------------------------
	// --------------------------operations----------------------------------------------------
	// out += lamda in
	friend void time(const double& lamda, const QMat& in, QMat& out);
	friend void timeLSign(const double& lamda, const QMat& in, QMat& out);

	// out += O * in
	friend void leftTime(const QMat& O, const QMat& in, QMat& out);

	// out += lamda O * in
	friend void leftTime(const double& lamda, const QMat& O, const QMat& in, QMat& out);
	friend void leftTimeLSign(const double& lamda, const QMat& O, const QMat& in, QMat& out);

	// out += in * O.t
	friend void rightTime(const QMat& O, const QMat& in, QMat& out);
	friend void rightTimeLSign(const double& lamda, const QMat& O, const QMat& in, QMat& out);

	// out += lamda in * O.t
	friend void rightTime(const double& lamda, const QMat& O, const QMat& in, QMat& out);


	// out += leftO * in * rightO
	friend void lrTime(const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out);

	// out += lamda leftO * in * rightO.t
	friend void lrTime(const double& lamda, const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out);
	friend void lrTimeLSign(const double& lamda, const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out);
	// --------------------------operations----------------------------------------------------
	//-----------------------------------------------------------------------------------------
	
	void zeros();
	int v2QMat(const double *f);
	int QMat2v(double *f)const;
	void print()const;
	void print(std::string)const;
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
