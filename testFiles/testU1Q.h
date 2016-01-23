#include <iostream>
#include "U1Q.h"


void testU1Q(){
	std::cout << "test U1Q\n";
	int q[] = {100,20};
	DMRGCat::U1Q q1(q);
	std::cout << q1 << std::endl;

	std::vector<int> qnos = {1,20};
	DMRGCat::U1Q q2(qnos);
	std::cout << q2 << "\n";
	std::cout << "\n\n";
	system("pause");
}