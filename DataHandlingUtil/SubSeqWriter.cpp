#include "SubSeqWriter.hpp"

class BaseSubSeqWriter : public ISubSeqWriter
{
public:
	void writeData(std::string& filename, std::vector<std::string> * subsequences)
	{
		std::ofstream outfile(filename.c_str());
		for (int idx = 0; idx < subsequences->size(); idx++)
		{
			outfile << (*subsequences)[idx] << std::endl;
		}
	}
};