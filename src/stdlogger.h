#include "logging.h"
#include <iostream>
 
 
using namespace std;
 
class StdLogger : public Logging{
public:
    virtual void log_info(string_view entry) override {
        cout << "[INFO] " << entry << endl;
    }
    virtual void log_warning(string_view entry) override {
        cout << "[WARNING] " << entry << endl;
    }
    virtual void log_error(string_view entry) override {
        cout << "[ERROR] " << entry << endl;
    }
};