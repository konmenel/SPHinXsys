#include <iostream>
#include <fstream>
#include <string>

/**
 * @brief Class to logs output to a file and stdout at the same time
 * 
 */
class LogOutput
{
public:
	LogOutput(const std::string &file_name)
	{
		log_file_.open(file_name);
		is_opened_ = true;
	}

	~LogOutput()
	{
		if (is_opened_) log_file_.close();
	}

	void close()
	{
		if (is_opened_)
		{
			log_file_.close();
			is_opened_ = false;
		}
	}

	// Overload the << operator to write to both stdout and log file
	template<typename T>
	LogOutput& operator<<(const T& t)
	{
		std::cout << t;
		log_file_ << t;

		return *this;
	}

	// Overload the << operator to handle std::endl
	LogOutput& operator<<(std::ostream& (*f)(std::ostream&))
	{
		f(std::cout);
		f(log_file_);

		return *this;
	}
	
protected:
	bool is_opened_;
	std::ofstream log_file_;
};
