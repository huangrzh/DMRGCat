	
#ifndef DMRGCATDATATYPE
#define DMRGCATDATATYPE

#include <unordered_map>



namespace DMRGCat{


	//hash function: integers smaller than 2^16 are safe;
	struct hash_func {
		int operator()(const std::pair < int, int > &addr) const{
			return ((addr.first << 16) ^ addr.second);
		}
	};
	

	//IntPair2IntHashMap									
	typedef std::unordered_map < std::pair < int, int >, int, hash_func> IntPair2IntHashMap;

	
	
	//IntPair
	typedef std::pair<int,int> IntPair;

	
}


#endif

