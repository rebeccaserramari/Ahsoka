#include "logging.h"
#include <iostream>
#include <fstream> 
 
using namespace std;
 
class FileLogger : public Logging{
public:
	virtual void log_info(string_view entry) override {
    	ofstream outfile;
    	outfile.open("logfile.log", ios_base::app);
      outfile << "[INFO] " << entry << endl;
    }
    virtual void log_warning(string_view entry) override {
        ofstream outfile;
    	outfile.open("logfile.log", ios_base::app);
      outfile << "[WARNING] " << entry << endl;
    }
    virtual void log_error(string_view entry) override {
        ofstream outfile;
    	outfile.open("logfile.log", ios_base::app);
      outfile << "[ERROR] " << entry << endl;
    }
};