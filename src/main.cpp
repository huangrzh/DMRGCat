#include "testFiles/testU1Q.h"
#include "testFiles/testBlockQBase.h"
#include "U1Q.h"


int main(){	
	int id00 = DMRGCat::getID({0,0});
	std::cout << id00 << std::endl;
	int id10 = DMRGCat::getID({ 1, 0 });
	std::cout << id10 << std::endl;
	int id11 = DMRGCat::getID({ 1, 1 });
	std::cout << id11 << std::endl;

	std::vector<double> coe = {1.0,-1.0};
	for (auto& x : coe){
		std::cout << x << std::endl;
	}
	//testU1Q();
	//testBlockQBase();
#ifdef VISUAL
	system("pause");
#endif
	return 0;
}