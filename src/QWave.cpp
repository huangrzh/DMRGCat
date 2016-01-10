#include "setting.h"
#include "QWave.h"



DMRGCat::QWave::QWave(){};
DMRGCat::QWave::~QWave(){};

DMRGCat::QWave::QWave(int totq, Block& sys, Block& m, Block& n, Block& env){
	setWave(totq,sys,m,n,env);
}


int DMRGCat::QWave::setWave(int totq, Block& sys, Block& m, Block& n, Block& env){
	Dim = 0;
	TotQNo = totq;

	if (MNQID2SysEnvNo.size() > 0){
		MNQID2SysEnvNo.clear();
		SysEnvQMat.clear();
	}

	int sysEnvNo = 0;
	//super block base and initial wave(whose elements all all zero);	
	for (const auto& xm : m.QSpace.SubQIDDim){
		for (const auto& xn : n.QSpace.SubQIDDim){
			int qmn = DMRGCat::getAddID(xm.first, xn.first);
			std::vector<std::pair<int, int>> lrdims;
			std::vector<std::pair<int, int>> lrqids;

			for (const auto& xs : sys.QSpace.SubQIDDim){
				for (const auto& xe : env.QSpace.SubQIDDim){
					int qse = DMRGCat::getAddID(xs.first, xe.first);
					DMRGCat::U1Q Q(DMRGCat::getAddID(qmn, qse));
					if (Q.getChargeNo() == TotQNo){
						Dim += (xs.second)*(xe.second);
						lrdims.push_back({ xs.second, xe.second });
						lrqids.push_back({ xs.first, xe.first });
					}
				}
			}

			if (lrdims.size() > 0){
				MNQID2SysEnvNo[{xm.first, xn.first }] = sysEnvNo;
				sysEnvNo++;
				DMRGCat::QMat qmat(lrqids, lrdims);
				SysEnvQMat.push_back(qmat);
			}
		}
	}
	return Dim;
}




void DMRGCat::QWave::oneBody(int flag, const QMat& O, QWave& out)const{
	int mnNo = MNQID2SysEnvNo.size();
	switch (flag){

		case 1:{
			for (const auto& x : MNQID2SysEnvNo){
				DMRGCat::leftTime(O, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
			break; 
		}


		case 2:{
			for (const auto& x : MNQID2SysEnvNo){
				auto xx = O.RQID2MatNo.find(x.first.first);
				if (xx != O.RQID2MatNo.end()){
					int realMNi = O.R2LID.at(xx->second);
					DMRGCat::time(O.SubMat.at(xx->second)(0,0), SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
				}				
			}
			break;
		}


		case 3:{
			for (const auto& x : MNQID2SysEnvNo){
				auto xx = O.RQID2MatNo.find(x.first.second);
				if (xx != O.RQID2MatNo.end()){
					int realMNi = O.R2LID.at(xx->second);
					DMRGCat::time(O.SubMat.at(xx->second)(0, 0), SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
				}
			}
			break;
		}


		case 4:{
			for (const auto& x : MNQID2SysEnvNo){
				DMRGCat::rightTime(O, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
			break;
		}

		default:
			break;
	}
	
}



void DMRGCat::QWave::oneBody(int flag, const double& lamda, const QMat& O, QWave& out)const{
	int mnNo = MNQID2SysEnvNo.size();
	
	switch (flag){
	case 1:{
		for (const auto& x : MNQID2SysEnvNo){
			DMRGCat::leftTime(lamda, O, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
		}
		break;
	}


	case 2:{
		for (const auto& x : MNQID2SysEnvNo){
			auto xx = O.RQID2MatNo.find(x.first.first);
			if (xx != O.RQID2MatNo.end()){
				int realMNi = O.R2LID.at(xx->second);
				DMRGCat::time(lamda*(O.SubMat.at(xx->second)(0, 0)), SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
		}
		break;
	}


	case 3:{
		for (const auto& x : MNQID2SysEnvNo){
			auto xx = O.RQID2MatNo.find(x.first.second);
			if (xx != O.RQID2MatNo.end()){
				int realMNi = O.R2LID.at(xx->second);
				DMRGCat::time(lamda*(O.SubMat.at(xx->second)(0, 0)), SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
		}
		break;
	}


	case 4:{
		for (const auto& x : MNQID2SysEnvNo){
			DMRGCat::rightTime(lamda, O, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
		}
		break;
	}


	default:
		break;
	}

}