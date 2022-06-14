#include <memory>
#include <string_view>
 
using namespace std;

//idea after the example from https://developer-blog.net/professionelles-loggen-unter-c/
class Logging {
public:
    virtual ~Logging() = default;
    virtual void log_info(string_view entry) = 0;
    virtual void log_warning(string_view entry) = 0;
    virtual void log_error(string_view entry) = 0;
};
 
using Logger = shared_ptr<Logging>;