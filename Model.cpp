#include "Model.h"


DMRGCat::Model::~Model(){}
DMRGCat::Model::Model(const std::string& model_name) {
	if (model_name == "FermionHubbard"){
		for (int i = 0; i < 2; i++){
			for (int j = 0; j < 2; j++){
				int qs[] = { i, j };
				DMRGCat::U1Q Q(qs);
				SingleSiteSpaceQID.push_back(Q.getID());
			}
		}

	}
}


void DMRGCat::Model::getSpaceQID(std::vector<int>& qids)const{
	qids = SingleSiteSpaceQID;
}
