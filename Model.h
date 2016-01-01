#ifndef MODEL_H
#define	MODEL_H

#include <string>
#include <vector>
#include "setting.h"
#include "U1Q.h"
namespace DMRGCat{

class Model{
public:
	Model(const std::string& model_name);
	void getSpaceQID(std::vector<int>&)const;
	~Model();		

private:		
	std::vector<int> SingleSiteSpaceQID;
};

}
#endif
