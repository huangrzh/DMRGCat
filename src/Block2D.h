#ifndef BLOCKTWOD_H
#define BLOCKTWOD_H


#include "QMat.h"
#include "Block.h"
#include "MulChain.h"

namespace DMRGCat{
	class Block2D : protected Block{
	public:
		Block2D(){}
		~Block2D(){}
		Block2D(const Block2D& var);
		Block2D(const Parameter& para, int site);		
		

		void initial2D(const Parameter& para, int site);
		void update2D(const Parameter& para, const MulChain& chain, const Block2D& added, const Block2D& old);
		void operator=(const Block2D& var);
		void reNorm2D(const BlockQBase& UBase, const QMat& reNormU);


		void save2D(std::ofstream& savefile)const;
		void save2D(int s)const;
		void load2D(std::ifstream& loadfile);
		void load2D(int s);
		void print2D()const;
		void print2D(std::string)const;


		friend class QWave;
		friend class SuperBlock2D;
	private:
		int SiteNo;
		int NoOfSites;
		std::unordered_map<int, int> InSite2OperatorNo;
		std::unordered_map<int, int> OutSite2OperatorNo;
		std::vector<QMat> OutQOCu;
		std::vector<QMat> OutQOCd;
		std::vector<QMat> OutQOCuDag;
		std::vector<QMat> OutQOCdDag;
		std::vector<QMat> InQOCu;
		std::vector<QMat> InQOCd;
		std::vector<QMat> InQOCuDag;
		std::vector<QMat> InQOCdDag;

		void clear2D();
	};
}


#endif