#ifndef QWAVE
#define QWAVE

#include <vector>
#include "datatype.h"
#include "QMat.h"
#include "Block.h"
#include <unordered_map>

namespace DMRGCat{


class QWave{
	

	public:	
		QWave();
		~QWave();
		QWave(int totq, Block& sys, Block& m, Block& n, Block& env);
		int setWave(int totq, Block& sys, Block& m, Block& n, Block& env);
		void getTruncU(BlockQBase& UBase, QMat& truncU);
	private:
		int Dim;
		int TotQNo;
		IntPair2IntHashMap MNQID2SysEnvNo;
		std::vector<QMat> SysEnvQMat;

		//operation:
		void oneBody(int flag, const QMat& O, QWave& out)const;
		void oneBody(int flag, const double& lamda, const QMat& O, QWave& out)const;

		void twoBody(int flagl, const QMat& Ol, int flagr, const QMat& Or, QWave& out)const;
		void twoBody(int flagl, const QMat& Ol, int flagr, const QMat& Or, const double& lamda, QWave& out)const;
			void SysMAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const;
			void SysNAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const;
			void SysEnvAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const;
			void MNAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const;
			void MEnvAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const;
			void NEnvAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const;

		void thrBody(int flag1, const QMat& O1, int flag2, const QMat& O2, int flag3, const QMat& O3, QWave& out)const;
		void thrBody(int flag1, const QMat& O1, int flag2, const QMat& O2, int flag3, const QMat& O3, const double& lamda, QWave& out)const;
			void SysMNAct(double lamda, const QMat& O1, const QMat& O2, const QMat& O3, QWave& out)const;
			void SysMEnvAct(double lamda, const QMat& O1, const QMat& O2, const QMat& O3, QWave& out)const;
			void SysNEnvAct(double lamda, const QMat& O1, const QMat& O2, const QMat& O3, QWave& out)const;
			void MNEnvAct(double lamda, const QMat& O1, const QMat& O2, const QMat& O3, QWave& out)const;

		void fourBody(const QMat& O1, const QMat& O2, const QMat& O3, const QMat& O4, QWave& out)const;
};
	
}


#endif

