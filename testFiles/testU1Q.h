#include <iostream>
#include "U1Q.h"


void testU1Q(){
	int q[] = {100,20};
	DMRGCat::U1Q q1(q);
	std::cout << q1 << std::endl;
}