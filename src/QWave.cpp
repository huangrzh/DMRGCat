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



void DMRGCat::QWave::twoBody(int flagl, const QMat& Ol, int flagr, const QMat& Or, QWave& out)const{
	double lamda = 1.0;
	int flag = flagl * 10 + flagr;

	switch (flag){
	case 12:{SysMAct(lamda, Ol, Or, out);break;}
	case 13:{SysNAct(lamda, Ol, Or, out); break;}		
	case 14:{SysEnvAct(lamda, Ol, Or, out); break;}
	case 23:{MNAct(lamda, Ol, Or, out);	break;}
	case 24:{MEnvAct(lamda, Ol, Or, out); break;}
	case 34:{NEnvAct(lamda, Ol, Or, out); break;}
	default:
		break;
	}
}


void DMRGCat::QWave::SysMAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const{
	if (Ol.getIsFermion()){
		for (const auto& x : MNQID2SysEnvNo){
			auto xm = Or.RQID2MatNo.find(x.first.first);
			if (xm != Or.RQID2MatNo.end()){
				auto newno = out.MNQID2SysEnvNo.find({ Ol.R2LID.at(x.first.first), x.first.second });
				if (newno != out.MNQID2SysEnvNo.end()){
					lamda *= Or.SubMat.at(Or.RQID2MatNo.at(x.first.first))(0, 0);
					lamda *= -(double)DMRGCat::getFermionSign(x.first.first);
					DMRGCat::leftTime(lamda, Ol, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
				}
			}
		}
	}
	else{
		for (const auto& x : MNQID2SysEnvNo){
			auto xm = Or.RQID2MatNo.find(x.first.first);
			if (xm != Or.RQID2MatNo.end()){
				lamda = Or.SubMat.at(Or.RQID2MatNo.at(x.first.first))(0, 0);
				DMRGCat::leftTime(lamda, Ol, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
		}
	}
}


void DMRGCat::QWave::SysNAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const{
	if (Ol.getIsFermion()){
		for (const auto& x : MNQID2SysEnvNo){
			auto xn = Or.RQID2MatNo.find(x.first.second);			
			if (xn != Or.RQID2MatNo.end()){
				auto newno = out.MNQID2SysEnvNo.find({ x.first.first, Or.R2LID.at(xn->first)});
				if (newno != out.MNQID2SysEnvNo.end()){
					lamda *= Or.SubMat.at(Or.RQID2MatNo.at(x.first.first))(0, 0);
					DMRGCat::leftTimeLSign(lamda, Ol, SysEnvQMat.at(x.second), out.SysEnvQMat.at(newno->second));
				}
			}
		}
	}
	else{
		for (const auto& x : MNQID2SysEnvNo){
			auto xm = Or.RQID2MatNo.find(x.first.first);
			if (xm != Or.RQID2MatNo.end()){
				lamda *= Or.SubMat.at(Or.RQID2MatNo.at(x.first.first))(0, 0);
				DMRGCat::leftTime(lamda, Ol, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
		}
	}
}


void DMRGCat::QWave::SysEnvAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const{
	if (Ol.getIsFermion()){
		for (const auto& x : MNQID2SysEnvNo){
			lamda *= (double)DMRGCat::getFermionSign(x.first.second);
			DMRGCat::lrTimeLSign(lamda, Ol, Or, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
		}
	}
	else{
		for (const auto& x : MNQID2SysEnvNo){
			lamda *= (double)DMRGCat::getFermionSign(x.first.second);
			DMRGCat::lrTime(lamda, Ol, Or, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
		}
	}
}




void DMRGCat::QWave::MNAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const{
	if (Ol.getIsFermion()){
		for (const auto& x : MNQID2SysEnvNo){
			auto xn = Or.RQID2MatNo.find(x.first.second);
			auto xm = Ol.RQID2MatNo.find(x.first.first);
			
			if (xm != Ol.RQID2MatNo.end() && xn != Or.RQID2MatNo.end()){
				auto newno = out.MNQID2SysEnvNo.find({ Ol.R2LID.at(x.first.first), Or.R2LID.at(x.first.second)});
				if (newno != out.MNQID2SysEnvNo.end()){
					lamda *= Or.SubMat.at(Or.RQID2MatNo.at(x.first.second))(0, 0);
					lamda *= Ol.SubMat.at(Ol.RQID2MatNo.at(x.first.first))(0, 0);
					lamda *= (double)DMRGCat::getFermionSign(x.first.first);
					DMRGCat::timeLSign(lamda, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
				}
			}
		}
	}
	else{
		for (const auto& x : MNQID2SysEnvNo){
			auto xn = Or.RQID2MatNo.find(x.first.second);
			auto xm = Ol.RQID2MatNo.find(x.first.first);
			if (xm != Ol.RQID2MatNo.end() && xn!=Or.RQID2MatNo.end()){
				lamda *= Or.SubMat.at(Or.RQID2MatNo.at(x.first.second))(0, 0);
				lamda *= Ol.SubMat.at(Ol.RQID2MatNo.at(x.first.first))(0, 0);
				DMRGCat::time(lamda,  SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
		}
	}
}




void DMRGCat::QWave::MEnvAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const{
	if (Ol.getIsFermion()){
		for (const auto& x : MNQID2SysEnvNo){
			auto xm = Ol.RQID2MatNo.find(x.first.first);
			if (xm != Ol.RQID2MatNo.end()){
				auto newno = out.MNQID2SysEnvNo.find({ Ol.R2LID.at(x.first.first), x.first.second });
				if (newno != out.MNQID2SysEnvNo.end()){
					lamda *= Ol.SubMat.at(Ol.RQID2MatNo.at(x.first.first))(0, 0);
					lamda *= (double)DMRGCat::getFermionSign(x.first.first, x.first.second);
					DMRGCat::rightTimeLSign(lamda, Or, SysEnvQMat.at(x.second), out.SysEnvQMat.at(newno->second));
				}
			}
		}
	}
	else{
		for (const auto& x : MNQID2SysEnvNo){
			auto xm = Ol.RQID2MatNo.find(x.first.first);
			if (xm != Ol.RQID2MatNo.end()){
				lamda *= Ol.SubMat.at(Ol.RQID2MatNo.at(x.first.first))(0, 0);
				DMRGCat::rightTime(lamda, Or, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
		}
	}
}


void DMRGCat::QWave::NEnvAct(double lamda, const QMat& Ol, const QMat& Or, QWave& out)const{
	if (Ol.getIsFermion()){
		for (const auto& x : MNQID2SysEnvNo){
			auto xn = Ol.RQID2MatNo.find(x.first.second);
			if (xn != Ol.RQID2MatNo.end()){
				auto newno = out.MNQID2SysEnvNo.find({ x.first.first, Ol.R2LID.at(x.first.second) });
				if (newno != out.MNQID2SysEnvNo.end()){
					lamda *= Ol.SubMat.at(Ol.RQID2MatNo.at(x.first.second))(0, 0);
					lamda *= (double)DMRGCat::getFermionSign(x.first.first, x.first.second);
					DMRGCat::rightTimeLSign(lamda, Or, SysEnvQMat.at(x.second), out.SysEnvQMat.at(newno->second));
				}
			}
		}
	}
	else{
		for (const auto& x : MNQID2SysEnvNo){
			auto xn = Ol.RQID2MatNo.find(x.first.second);
			if (xn != Ol.RQID2MatNo.end()){
				lamda *= Ol.SubMat.at(Ol.RQID2MatNo.at(x.first.second))(0, 0);
				DMRGCat::rightTime(lamda, Or, SysEnvQMat.at(x.second), out.SysEnvQMat.at(x.second));
			}
		}
	}
}



void DMRGCat::QWave::thrBody(int flag1, const QMat& O1, int flag2, const QMat& O2, int flag3, const QMat& O3, QWave& out)const{
	double lamda = 1.0;
	int flag = flag1 * 100 + flag2 * 10 + flag3;
	
	switch (flag){
	case 123:{SysMNAct(lamda, O1, O2, O3, out); break; }
	case 124:{SysMEnvAct(lamda, O1, O2, O3, out); break; }
	case 234:{MNEnvAct(lamda, O1, O2, O3, out); break; }
	case 134:{SysNEnvAct(lamda, O1, O2, O3, out); break; }
	default:
		break;
	}
	
}


void DMRGCat::QWave::thrBody(int flag1, const QMat& O1, int flag2, const QMat& O2, int flag3, const QMat& O3, const double& lamda, QWave& out)const{

}


void DMRGCat::QWave::fourBody(const QMat& O1, const QMat& O2, const QMat& O3, const QMat& O4, QWave& out)const{

}