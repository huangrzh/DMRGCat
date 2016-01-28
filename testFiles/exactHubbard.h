#include <armadillo>
#include <bitset>
#include <map>









//single site : |up, down>
//|i_L, i_L-1, i_L-2,............>
// |0,0> - 1
// |0,1> - 2
// |1,0> - 3
// |1,1> - 4
void exactHubbard(){
	std::cout << __FUNCTION__ << std::endl;

	const int L = 2;
	const int Nu = 1;
	const int Nd = 1;

	struct pairBitSetCMP {
		bool operator() (const std::pair<std::bitset<L>, std::bitset<L>>& lhs, const std::pair<std::bitset<L>, std::bitset<L>>& rhs) const{
			bool var = false;
			if (lhs.second.to_ulong() > rhs.second.to_ulong()){
				var = true;
			}
			else{
				if (lhs.first.to_ulong() > rhs.first.to_ulong()){
					var = true;
				}
			}
			//std::cout << "var = " << var << "\n";
			return  var;
		}
	};

	
	std::map<std::pair<std::bitset<L>,std::bitset<L>>, int, pairBitSetCMP> Maps;
	
	int nos = 0;
	unsigned int maxs = L*L;
	for (unsigned int nu = 0; nu < maxs; nu++){
		for (unsigned int nd = 0; nd < maxs; nd++){
			std::bitset<L> su(nu);
			std::bitset<L> sd(nd);
			if (su.count() == L / 2 && sd.count() == L / 2){
				//std::cout << su << ", " << sd << ", " << nos << "\n";
				Maps[{su, sd}] = nos;
				nos++;
			}
		}
	}
	

	for (const auto& x : Maps){
		std::cout << x.first.first << ", " << x.first.second << ", " << x.second << "\n";
		//std::cout << Maps.at(x.first) << std::endl;
	}

	int size = Maps.size();
	std::cout << "size = " << size << "\n";
	arma::mat H(size,size);
	H.zeros();
	for (const auto& x : Maps){
		std::cout << x.first.first << ", " << x.first.second << ", " << x.second << "\n";
		//up
		for (unsigned int ith = 1; ith < L; ith++){
			std::cout << "up ith = " << ith << "\n";
			if (x.first.first[L - ith-1] == 1 && x.first.first[L - ith] == 0){			
				std::bitset<L> lu = x.first.first;
				lu[L - ith - 1] = 0;
				lu[L - ith] = 1;
				int rowth = Maps.at({ lu, x.first.second });
				if (x.first.second[L - ith - 1] == 1){
					H(rowth, x.second) += -1.0;
				}
				else{
					H(rowth, x.second) += 1.0;
				}
			}			
		}
		//down
		for (unsigned int ith = 1; ith < L; ith++){
			std::cout << "down ith = " << ith << "\n";
			if (x.first.second[L - ith - 1] == 1 && x.first.second[L - ith] == 0){
				std::bitset<L> ld = x.first.second;
				ld[L - ith - 1] = 0;
				ld[L - ith] = 1;
				std::cout << x.first.first << ", " << x.first.second << ", " << ld << "\n";
				int rowth = Maps.at({ x.first.first, ld });	
				std::cout << "here\n";
				system("pause");
				if (x.first.first[L - ith - 1] == 1){
					H(rowth, x.second) += -1.0;
				}
				else{
					H(rowth, x.second) += 1.0;
				}
			}
		}
	}
	
	H = H + H.t();
	H.print("H");
	system("pause");
	arma::vec eigs = arma::eig_sym(H);
	eigs.print("eigs");
}