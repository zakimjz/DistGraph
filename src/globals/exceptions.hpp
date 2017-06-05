#ifndef __EXCEPTIONS_HPP__
#define __EXCEPTIONS_HPP__

#include <stdexcept>

class method_unimplemented : public std::runtime_error {
public:
  method_unimplemented(const char *err) : runtime_error(std::string("Method ") + err + "not implemented") {
  }
};

#endif

