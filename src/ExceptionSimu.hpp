#ifndef EXCEPTION_SIMU_HPP
#define EXCEPTION_SIMU_HPP

#include <iostream>
#include <exception>

class ExceptionSimu : public std::exception {
private : 
  virtual const char* what() const throw()
  {
    return "Exception Simulation";
  }
};

#endif
