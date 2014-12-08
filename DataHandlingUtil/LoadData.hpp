#ifndef LOAD_DATA_HPP_
#define LOAD_DATA_HPP_

// standard includes
#include <fstream>
#include <string>
#include <stdexcept>

class IDataLoader
{
	public:
		virtual std::string loadData(std::string& filename) = 0;
};

#endif /* LOAD_DATA_HPP_ */