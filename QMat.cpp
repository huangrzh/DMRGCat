#include "QMat.h"


DMRGCat::QMat::QMat(){}
DMRGCat::QMat::~QMat(){}






DMRGCat::QMat::QMat(const std::vector<std::pair<int, int>>& LRIDs, const std::vector<double>& coe){
	try{
		if (LRIDs.size() != coe.size()){
			throw std::runtime_error("Error in QMat LRIDs.size!=coe.size");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
	int submatno = 0;
	for (int i = 0; i < LRIDs.size(); i++){
		arma::mat temp_mat(1,1);
		temp_mat = coe.at(i);
		SubMat.push_back(temp_mat);
		RQID2MatNo[LRIDs[i].second] = submatno;
		LQID2MatNo[LRIDs[i].first] = submatno;
		submatno++;
	}
}



DMRGCat::QMat::QMat(const QMat& qmat1, const QMat& qmat2, const BlockQBase& bigBase){
	int no = 0;
	for (const auto& x : qmat1.R2LID){
		for (const auto& y : qmat2.R2LID){
			int newRID = DMRGCat::getAddID(x.first,y.first);
			int newLID = DMRGCat::getAddID(x.second,y.second);
			arma::mat tempMat1 = arma::kron(qmat1.SubMat.at(qmat1.RQID2MatNo(x.first)), qmat2.SubMat.at(qmat2.RQID2MatNo.at(y.first)));
			if (R2LID.find(newRID) == R2LID.end()){
				R2LID[newRID] = newRID;
				arma::mat tempMat(bigBase.SubQIDDim.at(newLID), bigBase.SubQIDDim.at(newRID));				
				matCompress(tempMat, bigBase.StartDim.at(std::pair<int,int>(x.second,y.second)), tempMat1, bigBase.StartDim.at(std::pair<int,int>(x.first,y.first)));
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

	for (const auto& x : QMatvar.R2LID){
		output << "Q1:" << DMRGCat::U1Q(x.second) << "Q2:" << DMRGCat::U1Q(x.second);
		int no = QMatvar.RQID2MatNo.at(x.first);
		output << "Matrix: " << QMatvar.SubMat.at(no).n_rows << " x " << QMatvar.SubMat.at(no).n_cols << std::endl;

		//if (QMatvar.SubMat.at(no).n_rows <= 9 && QMatvar.SubMat.at(no).n_cols <= 9)
		//output << QMatvar.SubMat.at(no);
	}

	output << std::endl << std::endl;

	return output;
}


void DMRGCat::QMat::clear(){
	R2LID.clear();
	RQID2MatNo.clear();
	LQID2MatNo.clear();
	SubMat.clear();
}

void DMRGCat::QMat::save(std::ofstream& savefile) const{
	size_t NumOfBlock = SubMat.size();
	savefile.write((const char*)&NumOfBlock, sizeof(int));
	for (const auto& x : R2LID){
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
	arma::mat matrix;
	for (int i = 0; i < NumOfBlock; i++){
		loadfile.read((char*)&(x), sizeof(int));
		loadfile.read((char*)&(y), sizeof(int));
		R2LID[x] = y;
	}

	for (int i = 0; i < NumOfBlock; i++){
		loadfile.read((char*)&(x), sizeof(int));
		loadfile.read((char*)&(y), sizeof(int));
		RQID2MatNo[x] = y;
		LQID2MatNo[R2LID.at(x)] = y;
	}

	for (int i = 0; i < NumOfBlock; i++){
		loadfile.read((char*)&(x), sizeof(int));
		loadfile.read((char*)&(y), sizeof(int));
		R2LID[x] = y;
	}
		
	for (int i = 0; i < NumOfBlock; i++){
		try{
			bool loadok = matrix.load(loadfile, arma::arma_binary);
			if (!loadok){
				throw std::runtime_error("Error in QMat::load");
			}
		}
		catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
		}
	}
}