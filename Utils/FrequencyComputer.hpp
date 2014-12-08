#ifndef FREQ_COMP_HPP_
#define FREQ_COMP_HPP_

#include <map>
#include <stdexcept>

class IFrequencyComputer
{
	public:
		virtual std::string getSymbols() = 0;
		virtual double * computeFrequency(std::string& sequence) = 0;
		virtual void initChrCounts(std::map <char, long>& chrCounts) = 0;
		virtual void countCharsInSeq(std::string& sequence, std::map <char, long>& chrCounts) = 0;
		virtual void increaseCountChar(char chrInSeq, std::map <char, long>& chrCounts) = 0;
		virtual void decreaseCountChar(char chrInSeq, std::map <char, long>& chrCounts) = 0;
		virtual double * determineFrequencies(std::map <char, long>& chrCounts, int sequenceLength) = 0;
};

class FrequencyComputer : public IFrequencyComputer
{
	public:
		FrequencyComputer(std::string symbols) : symbols(symbols)
		{
		}
		virtual std::string getSymbols()
		{
			return symbols;
		}
		virtual double * computeFrequency(std::string&);
		virtual void initChrCounts(std::map <char, long>&);
		virtual void countCharsInSeq(std::string&, std::map <char, long>&);
		virtual void increaseCountChar(char, std::map <char, long>&);
		virtual void decreaseCountChar(char, std::map <char, long>&);
		virtual double * determineFrequencies(std::map <char, long>&, int);
		
	private:
		void throwErrorMsg(char currentChr)
		{
			std::string errmsg("The input sequence contains such symbol '");
			errmsg.append(std::string(1, currentChr));
			errmsg.append("' which is not valid.\n");
			throw std::runtime_error(errmsg);
		}

		void determineFrequencies(double * frequencies, std::map <char, long>& chrCounts, int sequenceLength)
		{
			for (int idx = 0; idx < symbols.length(); idx++)
			{
				frequencies[idx] = (double)chrCounts[symbols[idx]] / (double)sequenceLength;
			}
		}

		std::string symbols;
};

#endif /* FREQ_COMP_HPP_ */