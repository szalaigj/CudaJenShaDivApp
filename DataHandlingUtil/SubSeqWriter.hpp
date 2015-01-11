#ifndef SUB_SEQ_WRITER_HPP_
#define SUB_SEQ_WRITER_HPP_

// standard includes
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

class ISubSeqWriter
{
public:
	virtual void writeData(std::string& filename, std::vector<std::string> * subsequences) = 0;
};

#endif /* SUB_SEQ_WRITER_HPP_ */