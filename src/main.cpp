#include <ostream>
#include <istream>
#include "../testFiles/testU1Q.h"
#include "../testFiles/testBlockQBase.h"
#include "../testFiles/testQMat.h"
#include "../testFiles/testBlock.h"
#include "../testFiles/testQWave.h"
#include "U1Q.h"
#include "Parameter.h"
#include "Block.h"


int main(){	
	//testU1Q();
	//testBlockQBase();	
	//testQMat();
	//testBlock();
	testQWave();


#ifdef VISUAL
	system("pause");
#endif
	return 0;
}