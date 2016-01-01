#include "BlockQBase.h"



DMRGCat::BlockQBase::BlockQBase(){}
DMRGCat::BlockQBase::~BlockQBase(){}



DMRGCat::BlockQBase::BlockQBase(const DMRGCat::Model& model){
	std::vector<int> qids;
	model.getSpaceQID(qids);
	for (const auto& x : qids){
		SubQIDDim[x] = 1;
	}
}



//base b1 and b2 form bigger base *this;
//In DMRG, b1 for block left or right, b2 for added single site
//				  %------------------%	
//                | NewSys = Sys + m |
//				  %------------------%	
DMRGCat::BlockQBase::BlockQBase(const DMRGCat::BlockQBase &b1, const DMRGCat::BlockQBase &b2){
	for (const auto& x1 : b1.SubQIDDim){		
		for (const auto& x2 : b2.SubQIDDim){
			int id = DMRGCat::getAddID(x1.first,x2.first);
			auto find_id = SubQIDDim.find(id);
			if (find_id != SubQIDDim.end()){
				StartDim[std::pair<int, int>(x1.first, x2.first)] = find_id->second;
				find_id->second += x1.second * x2.second;
			}
			else{
				StartDim[std::pair<int, int>(x1.first, x2.first)] = 0;
				SubQIDDim[id] = x1.second * x2.second;				
			}
		}
	}
}



void DMRGCat::getKronOrder(const DMRGCat::BlockQBase &b1, const DMRGCat::BlockQBase &b2, std::map<std::pair<int, int>, int>& start_dim){
	start_dim.clear();
	std::unordered_map<int, int> qid_dim;
	for (const auto& x1 : b1.SubQIDDim){
		for (const auto& x2 : b2.SubQIDDim){
			int id = DMRGCat::getAddID(x1.first, x2.first);
			auto find_id = qid_dim.find(id);
			if (find_id != qid_dim.end()){
				start_dim[std::pair<int, int>(x1.first, x2.first)] = find_id->second;
				find_id->second += x1.second * x2.second;
			}
			else{
				start_dim[std::pair<int, int>(x1.first, x2.first)] = 0;
				qid_dim[id] = x1.second * x2.second;
			}
		}
	}
}



//Use truncating U, which is got from density matrix or wave function, to truncate space;
void DMRGCat::BlockQBase::truncate(DMRGCat::BlockQBase &ubase){
	std::vector<int> eraseqs;
	for (auto& x : SubQIDDim){
		auto itfind = ubase.SubQIDDim.find(x.first);
		if (itfind == ubase.SubQIDDim.end()){
			//It's important not to erase element insiede the for loop here;
			//Erase may(maybe not, but for sure we do it outside) cause rehash
			//of SubQIDDim which makes the for loop not right any longer;
			eraseqs.push_back(x.first);
		}
		else{
			x.second = itfind->second;
		}
	}

	for (const auto& x : eraseqs){
		SubQIDDim.erase(x);
	}
}



std::ostream& DMRGCat::operator<<(std::ostream& output, const DMRGCat::BlockQBase& BaseVar){
	for (const auto& x : BaseVar.SubQIDDim){
		DMRGCat::U1Q Q(x.first);
		output << Q << " -> " << x.second << std::endl;
	}
	return output;
}