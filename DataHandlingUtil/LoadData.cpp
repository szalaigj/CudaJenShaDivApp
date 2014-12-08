#include "LoadData.hpp"

class BaseDataLoader : public IDataLoader
{
	public:
		virtual std::string loadData(std::string& filename)
		{
			std::string result("");
			std::string line;
			std::ifstream inputfile(filename.c_str());
			if (inputfile.good())
			{
				while (std::getline(inputfile, line))
				{
					result += line;
				}
				inputfile.close();
			}
			else
			{
				std::string errmsg("File ");
				errmsg.append(filename);
				errmsg.append(" does not exist.\n");
				throw std::runtime_error(errmsg);
			}
			return result;
		}
};