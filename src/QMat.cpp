#include "U1Q.h"
#include "QMat.h"
#include "setting.h"

DMRGCat::QMat::QMat(){}
DMRGCat::QMat::~QMat(){}



DMRGCat::QMat::QMat(const QMat& var){
	IsFermion = var.IsFermion;
	RQID2MatNo = var.RQID2MatNo;
	LQID2MatNo = var.LQID2MatNo;
	R2LID = var.R2LID;
	LRID = var.LRID;
	SubMat = var.SubMat;
}


void DMRGCat::QMat::operator=(const QMat& var){
	IsFermion = var.IsFermion;
	RQID2MatNo = var.RQID2MatNo;
	LQID2MatNo = var.LQID2MatNo;
	R2LID = var.R2LID;
	LRID = var.LRID;
	SubMat = var.SubMat;
}


DMRGCat::QMat::QMat(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<std::pair<int, int>>& LRDims){
	IsFermion = DMRGCat::hasSign(LRIDs.at(0).first, LRIDs.at(0).second);
	zero(LRIDs,LRDims);
}


void DMRGCat::QMat::zero(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<std::pair<int, int>>& LRDims){
	clear();

	IsFermion = DMRGCat::hasSign(LRIDs.at(0).first, LRIDs.at(0).second);
	try{
		if (LRIDs.size() != LRDims.size()){
			throw std::runtime_error("Error in QMat LRIDs.size!=LRDims.size");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
	int submatno = 0;
	LRID = LRIDs;
	for (int i = 0; i < LRIDs.size(); i++){
		R2LID[LRIDs[i].second] = LRIDs[i].first;
		arma::mat temp_mat;
		temp_mat.zeros(LRDims.at(i).first, LRDims.at(i).second);
		SubMat.push_back(temp_mat);
		RQID2MatNo[LRIDs[i].second] = submatno;
		LQID2MatNo[LRIDs[i].first] = submatno;
		submatno++;
	}
}


DMRGCat::QMat::QMat(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<double>& coe){
	IsFermion = DMRGCat::hasSign(LRIDs.at(0).first, LRIDs.at(0).second);
	try{
		if (LRIDs.size() != coe.size()){
			throw std::runtime_error("Error in QMat LRIDs.size!=coe.size");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
	int submatno = 0;
	LRID = LRIDs;
	for (int i = 0; i < LRIDs.size(); i++){
		R2LID[LRIDs[i].second] = LRIDs[i].first;
		arma::mat temp_mat(1,1);
		temp_mat = coe.at(i);
		SubMat.push_back(temp_mat);
		RQID2MatNo[LRIDs[i].second] = submatno;
		LQID2MatNo[LRIDs[i].first] = submatno;
		submatno++;
	}
}



DMRGCat::QMat::QMat(const std::vector<std::pair<int, int>>& LRIDs){
	IsFermion = DMRGCat::hasSign(LRIDs.at(0).first, LRIDs.at(0).second);
	LRID = LRIDs;
	int submatno = 0;
	for (int i = 0; i < LRIDs.size(); i++){
		R2LID[LRIDs[i].second] = LRIDs[i].first;
		arma::mat temp_mat;
		temp_mat.eye(1, 1);
		SubMat.push_back(temp_mat);
		RQID2MatNo[LRIDs[i].second] = submatno;
		LQID2MatNo[LRIDs[i].first] = submatno;
		submatno++;
	}
}


DMRGCat::QMat::QMat(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase){
	IsFermion = (qmat1.IsFermion || qmat2.IsFermion) && (!qmat1.IsFermion || !qmat2.IsFermion);
	kron(qmat1, qmat2, bigBase);
}

//qmat1 x qmat2  ->  |qmat1.i, qmat2.i>
void DMRGCat::QMat::kron(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase){
	clear();
	IsFermion = (qmat1.IsFermion || qmat2.IsFermion) && (!qmat1.IsFermion || !qmat2.IsFermion);
	int no = 0;
#ifdef FERMION
	if (qmat2.getIsFermion()){
		for (const auto& x : qmat1.R2LID){
			for (const auto& y : qmat2.R2LID){
				int newRID = DMRGCat::getAddID(x.first, y.first);
				int newLID = DMRGCat::getAddID(x.second, y.second);
				arma::mat tempMat1;
				if (DMRGCat::getFermionSign(x.second) == -1){					
					tempMat1 = -arma::kron(qmat1.SubMat.at(qmat1.RQID2MatNo.at(x.first)), qmat2.SubMat.at(qmat2.RQID2MatNo.at(y.first)));
				}
				else{
					tempMat1 = arma::kron(qmat1.SubMat.at(qmat1.RQID2MatNo.at(x.first)), qmat2.SubMat.at(qmat2.RQID2MatNo.at(y.first)));
				}
				if (R2LID.find(newRID) == R2LID.end()){
					LRID.push_back({ newLID, newRID });
					R2LID[newRID] = newLID;
					arma::mat tempMat;
					tempMat.zeros(bigBase.SubQIDDim.at(newLID), bigBase.SubQIDDim.at(newRID));
					matCompress(tempMat, bigBase.StartDim.at(std::pair<int, int>(x.second, y.second)), tempMat1, bigBase.StartDim.at(std::pair<int, int>(x.first, y.first)));
					SubMat.push_back(tempMat);
					RQID2MatNo[newRID] = no;
					LQID2MatNo[newLID] = no;
					no++;
				}
				else{
					matCompress(SubMat.at(RQID2MatNo.at(newRID)), bigBase.StartDim.at(std::pair<int, int>(x.second, y.second)), tempMat1, bigBase.StartDim.at(std::pair<int, int>(x.first, y.first)));
				}
			}
		}
	}
	else{
		for (const auto& x : qmat1.R2LID){
			for (const auto& y : qmat2.R2LID){
				int newRID = DMRGCat::getAddID(x.first, y.first);
				int newLID = DMRGCat::getAddID(x.second, y.second);
				arma::mat tempMat1 = arma::kron(qmat1.SubMat.at(qmat1.RQID2MatNo.at(x.first)), qmat2.SubMat.at(qmat2.RQID2MatNo.at(y.first)));
				
				if (R2LID.find(newRID) == R2LID.end()){
					LRID.push_back({ newLID, newRID });
					R2LID[newRID] = newLID;
					arma::mat tempMat;
					tempMat.zeros(bigBase.SubQIDDim.at(newLID), bigBase.SubQIDDim.at(newRID));
					matCompress(tempMat, bigBase.StartDim.at(std::pair<int, int>(x.second, y.second)), tempMat1, bigBase.StartDim.at(std::pair<int, int>(x.first, y.first)));
					SubMat.push_back(tempMat);
					RQID2MatNo[newRID] = no;
					LQID2MatNo[newLID] = no;
					no++;
				}
				else{
					matCompress(SubMat.at(RQID2MatNo.at(newRID)), bigBase.StartDim.at(std::pair<int, int>(x.second, y.second)), tempMat1, bigBase.StartDim.at(std::pair<int, int>(x.first, y.first)));
				}
			}
		}
	}
#endif	
}


void DMRGCat::QMat::matCompress(arma::mat& newMat, int dimL, const arma::mat& old, int dimR){
	for (int a = 1; a<old.n_rows + 1; a++){
		for (int b = 1; b<old.n_cols + 1; b++){
			newMat(dimL + a - 1, dimR + b - 1) = old(a - 1, b - 1);
		}
	}
}



std::ostream& DMRGCat::operator<<(std::ostream& output, const DMRGCat::QMat& QMatvar){
	output << "----------------------------------" << std::endl;
	output << "NumOfBlock:\t" << QMatvar.SubMat.size() << std::endl;
	std::cout << "size = " << QMatvar.R2LID.size() << "\n";
	for (const auto& x : QMatvar.R2LID){
		output << "Q1:" << DMRGCat::U1Q(x.second) << "Q2:" << DMRGCat::U1Q(x.first);
		int no = QMatvar.RQID2MatNo.at(x.first);
		std::cout << "no = " << no << "\n";
		output << "Matrix: " << QMatvar.SubMat.at(no).n_rows << " x " << QMatvar.SubMat.at(no).n_cols << std::endl;

		if (QMatvar.SubMat.at(no).n_rows <= 9 && QMatvar.SubMat.at(no).n_cols <= 9){
			output << QMatvar.SubMat.at(no) << "\n";
		}
	}

	output << std::endl << std::endl;

	return output;
}



void DMRGCat::QMat::eyeQMat(const DMRGCat::BlockQBase& base){
	IsFermion = false;
	clear();
	int submatno = 0;
	for (const auto& x : base.SubQIDDim){
		LRID.push_back({ x.first, x.first });
		R2LID[x.first] = x.first;
		arma::mat temp_mat;		
		SubMat.push_back(temp_mat);
		SubMat.at(submatno).eye(x.second, x.second);
		RQID2MatNo[x.first] = submatno;
		LQID2MatNo[x.first] = submatno;
		submatno++;
	}
}




void DMRGCat::QMat::clear(){
	LRID.clear();
	R2LID.clear();
	RQID2MatNo.clear();
	LQID2MatNo.clear();
	SubMat.clear();
}

void DMRGCat::QMat::save(std::ofstream& savefile) const{
	unsigned int NumOfBlock = SubMat.size();
	savefile.write((const char*)&NumOfBlock, sizeof(int));
	for (const auto& x : LRID){
		savefile.write((const char*)&x.first, sizeof(int));
		savefile.write((const char*)&x.second, sizeof(int));
	}

	for (const auto& x : RQID2MatNo){
		savefile.write((const char*)&x.first, sizeof(int));
		savefile.write((const char*)&x.second, sizeof(int));
	}
	for (const auto& x : SubMat){
		try{ 
			bool saveok = x.save(savefile, arma::arma_binary); 
			if (!saveok){
				throw std::runtime_error("Error in QMat::save");
			}
		}
		catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
		}
	}
}


void DMRGCat::QMat::load(std::ifstream& loadfile){
	clear();
	int NumOfBlock;
	loadfile.read((char*)&NumOfBlock, sizeof(int));	
	
	int x,y;
	for (int i = 0; i < NumOfBlock; i++){
		loadfile.read((char*)&(x), sizeof(int));
		loadfile.read((char*)&(y), sizeof(int));
		R2LID[y] = x;
		LRID.push_back({ x, y });
	}

	IsFermion = DMRGCat::hasSign(LRID.at(0).first, LRID.at(0).second);

	for (int i = 0; i < NumOfBlock; i++){
		loadfile.read((char*)&(x), sizeof(int));
		loadfile.read((char*)&(y), sizeof(int));
		RQID2MatNo[x] = y;
		LQID2MatNo[R2LID.at(x)] = y;
	}

	
	arma::mat matrix;
	SubMat = std::vector<arma::mat>(NumOfBlock,matrix);
	for (int i = 0; i < NumOfBlock; i++){
		try{
			bool loadok = SubMat.at(i).load(loadfile, arma::arma_binary);
			if (!loadok){
				throw std::runtime_error("Error in QMat::load");
			}
		}
		catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
		}
	}
	
}



void DMRGCat::QMat::trunc(const BlockQBase& UBase, const QMat& truncU){
	std::unordered_map<int, int> saveR2No = RQID2MatNo;
	std::vector<std::pair<int, int>> saveLRID = LRID;
	std::vector<arma::mat> saveMat = SubMat;
	std::map<std::pair<int, int>, int> truncLRID;
	clear();


	for (const auto& x : saveLRID){
		if (UBase.SubQIDDim.find(x.first) == UBase.SubQIDDim.end() || UBase.SubQIDDim.find(x.second) == UBase.SubQIDDim.end()){
			truncLRID[x]=1;
		} 
	}

	int subno = 0;
	for (const auto&x : saveLRID){
		if (truncLRID.find(x) == truncLRID.end()){
			int no = saveR2No.at(x.second);
			SubMat.push_back(saveMat.at(no));
			RQID2MatNo[x.second] = subno;
			LQID2MatNo[x.first] = subno;
			R2LID[x.second] = x.first;
			LRID.push_back(x);
			subno++;
		}
	}

	int no = 0;
	int rno = 0;
	int lno = 0;
	for (const auto&x : RQID2MatNo){
		no = x.second;
		lno = truncU.RQID2MatNo.at(R2LID.at(x.first));
		rno = truncU.RQID2MatNo.at(x.first);
		SubMat.at(no) = (truncU.SubMat.at(lno).t()) * (SubMat.at(no)) * (truncU.SubMat.at(rno));
	}
}


void DMRGCat::QMat::trans(const QMat& m){
	clear();
#ifdef FERMION
	IsFermion = m.IsFermion;
#endif
	int size = m.LRID.size();
	for (int i = 0; i < size; i++){
		LRID.push_back({m.LRID.at(i).second, m.LRID.at(i).first});
		R2LID[m.LRID.at(i).first] = m.LRID.at(i).second;
		SubMat.push_back(m.SubMat.at(i));
	}
	LQID2MatNo = m.RQID2MatNo;
	RQID2MatNo = m.LQID2MatNo;
}

void DMRGCat::QMat::trans(){
	std::unordered_map<int, int> r2l;
	std::vector<std::pair<int, int>> lr;
	int size = LRID.size();
	for (int i = 0; i < size; i++){
		lr.push_back({ LRID.at(i).second, LRID.at(i).first });
		r2l[LRID.at(i).first] = LRID.at(i).second;
		SubMat.at(i) = SubMat.at(i).t();
	}

	LRID = lr;
	R2LID = r2l;
	auto x = LQID2MatNo;
	LQID2MatNo = RQID2MatNo;
	RQID2MatNo = x;
}



void DMRGCat::QMat::add(const QMat& added, const BlockQBase& space){
	for (const auto& x : space.SubQIDDim){
		auto findx1 = RQID2MatNo.find(x.first);
		auto findx2 = added.RQID2MatNo.find(x.first);
		if (findx1 != RQID2MatNo.end() && findx2!=added.RQID2MatNo.end()){
			SubMat.at(findx1->second) += added.SubMat.at(findx2->second);
		}

		else if (findx1 == RQID2MatNo.end() && findx2 != added.RQID2MatNo.end()){
			int lqid = added.R2LID.at(x.first);
			LRID.push_back({ lqid, x.first });
			R2LID[x.first] = lqid;
			SubMat.push_back(added.SubMat.at(added.RQID2MatNo.at(x.first)));
			RQID2MatNo[x.first] = SubMat.size() - 1;
			LQID2MatNo[lqid] = SubMat.size() - 1;
		}
		else{ ; }
	}
}



bool DMRGCat::QMat::getIsFermion()const{
	return IsFermion;
}



void DMRGCat::QMat::time(double lamda){
	for (auto& x : SubMat){
		x *= lamda;
	}
}





//-----------------------------------------------------------------------------------
//------------------------------QMat operations--------------------------------------
//-----------------------------------------------------------------------------------

void DMRGCat::time(const double& lamda, const QMat& in, QMat& out){
	int matNo = in.SubMat.size();
	for (int i = 0; i < matNo; i++){
		out.SubMat.at(i) += lamda * in.SubMat.at(i);
	}
}





void DMRGCat::timeLSign(const double& lamda, const QMat& in, QMat& out){
	int matNo = in.SubMat.size();
	double coe = 1.0;
	for (int i = 0; i < matNo; i++){
		coe = lamda * (double)DMRGCat::getFermionSign(in.LRID.at(i).first);
		out.SubMat.at(i) += coe * in.SubMat.at(i);
	}
}




// out += O * in
void DMRGCat::leftTime(const QMat& O, const QMat& in, QMat& out){
	for (const auto& inx : in.LQID2MatNo){
		const auto x = O.RQID2MatNo.find(inx.first);
		if (x != O.RQID2MatNo.end()){
			out.SubMat.at(out.LQID2MatNo.at(O.R2LID.at(x->first))) += O.SubMat.at(x->second) * in.SubMat.at(inx.second);
		}
	}
}



void DMRGCat::leftTime(const double& lamda, const QMat& O, const QMat& in, QMat& out){
	for (const auto& inx : in.LQID2MatNo){
		const auto x = O.RQID2MatNo.find(inx.first);
		if (x != O.RQID2MatNo.end()){
			out.SubMat.at(out.LQID2MatNo.at(O.R2LID.at(x->first))) += lamda * O.SubMat.at(x->second) * in.SubMat.at(inx.second);
		}
	}
}


void DMRGCat::leftTimeLSign(const double& lamda, const QMat& O, const QMat& in, QMat& out){
	double coe = 1.0;
	for (const auto& inx : in.LQID2MatNo){
		const auto x = O.RQID2MatNo.find(inx.first);
		if (x != O.RQID2MatNo.end()){
			coe = lamda * (double)DMRGCat::getFermionSign(inx.first);
			out.SubMat.at(out.LQID2MatNo.at(O.R2LID.at(x->first))) += coe * O.SubMat.at(x->second) * in.SubMat.at(inx.second);
		}
	}
}




// out += in * O.t
void DMRGCat::rightTime(const QMat& O, const QMat& in, QMat& out){
	for (const auto& inx : in.RQID2MatNo){
		const auto x = O.RQID2MatNo.find(inx.first);
		if (x != O.RQID2MatNo.end()){
			out.SubMat.at(out.RQID2MatNo.at(O.R2LID.at(x->first))) += in.SubMat.at(inx.second) * O.SubMat.at(x->second).t();
		}
	}
}



void DMRGCat::rightTime(const double& lamda, const QMat& O, const QMat& in, QMat& out){
	for (const auto& inx : in.RQID2MatNo){
		const auto x = O.RQID2MatNo.find(inx.first);
		if (x != O.RQID2MatNo.end()){
			out.SubMat.at(out.RQID2MatNo.at(O.R2LID.at(x->first))) += lamda * in.SubMat.at(inx.second) * O.SubMat.at(x->second).t();
		}
	}
}


void DMRGCat::rightTimeLSign(const double& lamda, const QMat& O, const QMat& in, QMat& out){
	double coe = 1.0;
	for (const auto& inx : in.RQID2MatNo){
		const auto x = O.RQID2MatNo.find(inx.first);
		if (x != O.RQID2MatNo.end()){
			coe = lamda * (double)DMRGCat::getFermionSign(in.R2LID.at(inx.first));
			out.SubMat.at(out.RQID2MatNo.at(O.R2LID.at(x->first))) += coe * in.SubMat.at(inx.second) * O.SubMat.at(x->second).t();
		}
	}
}



// out += leftO * in * rightO
void DMRGCat::lrTime(const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out){
	for (const auto& inx : in.LRID){
		const auto lx = leftO.RQID2MatNo.find(inx.first);
		const auto rx = rightO.RQID2MatNo.find(inx.second);
		if (lx != leftO.RQID2MatNo.end() && rx != rightO.RQID2MatNo.end()){
			out.SubMat.at(out.RQID2MatNo.at(rightO.R2LID.at(rx->first))) += leftO.SubMat.at(lx->second) * in.LQID2MatNo.at(inx.first) * rightO.SubMat.at(rx->second).t();
		}
	}
}



void DMRGCat::lrTime(const double& lamda, const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out){
	for (const auto& inx : in.LRID){
		const auto lx = leftO.RQID2MatNo.find(inx.first);
		const auto rx = rightO.RQID2MatNo.find(inx.second);
		if (lx != leftO.RQID2MatNo.end() && rx != rightO.RQID2MatNo.end()){
			out.SubMat.at(out.RQID2MatNo.at(rightO.R2LID.at(rx->first))) += lamda * leftO.SubMat.at(lx->second) * in.LQID2MatNo.at(inx.first) * rightO.SubMat.at(rx->second).t();
		}
	}
}


void DMRGCat::lrTimeLSign(const double& lamda, const QMat& leftO, const QMat& rightO, const QMat& in, QMat& out){
	double coe = 1.0;
	for (const auto& inx : in.LRID){
		const auto lx = leftO.RQID2MatNo.find(inx.first);
		const auto rx = rightO.RQID2MatNo.find(inx.second);
		if (lx != leftO.RQID2MatNo.end() && rx != rightO.RQID2MatNo.end()){
			coe = lamda * (double)DMRGCat::getFermionSign(inx.first);
			out.SubMat.at(out.RQID2MatNo.at(rightO.R2LID.at(rx->first))) += lamda * leftO.SubMat.at(lx->second) * in.LQID2MatNo.at(inx.first) * rightO.SubMat.at(rx->second).t();
		}
	}
}
// --------------------------operations----------------------------------------------------
//-----------------------------------------------------------------------------------------


void DMRGCat::QMat::zeros(){
	for (auto& x : SubMat){
		x.zeros();
	}
}

int DMRGCat::QMat::v2QMat(const double* f){
	int id = 0;
	for (auto& x:SubMat){
		memcpy(x.memptr(), &f[id], x.n_elem*sizeof(double));
		id += x.n_elem;
	}
	return id;
}

int DMRGCat::QMat::QMat2v(double* f) const{
	int id = 0;
	for (const auto& x:SubMat)	{
		memcpy(&f[id], x.memptr(), x.n_elem*sizeof(double));
		id += x.n_elem;
	}
	return id;
}




void DMRGCat::QMat::print()const{
	std::cout << *this;
}



void DMRGCat::QMat::print(std::string s)const{
	std::cout << s << "\n";
	print();
}