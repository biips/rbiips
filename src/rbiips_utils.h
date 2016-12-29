#ifndef RBIIPS_UTILS_H_
#define RBIIPS_UTILS_H_

#include <Rcpp.h>
#include <vector>
#include <Console.hpp>
#include <common/Error.hpp>
#include "Rostream.h"
#include <common/NumArray.hpp>

#ifndef BEGIN_RBIIPS
#define BEGIN_RBIIPS BEGIN_RCPP
#endif

#ifndef CATCH_RBIIPS
#define CATCH_RBIIPS                                                    \
  }                                                                     \
  catch (Biips::RuntimeError & except)                                  \
  {                                                                     \
    forward_exception_to_r(except);                                          \
  }                                                                     \
  catch (Biips::LogicError & except)                                    \
  {                                                                     \
    forward_exception_to_r(except);
#endif

#ifndef VOID_END_RBIIPS
#define VOID_END_RBIIPS CATCH_RBIIPS VOID_END_RCPP
#endif

#ifndef END_RBIIPS
#define END_RBIIPS CATCH_RBIIPS END_RCPP
#endif

using namespace Biips;
using std::endl;

extern Size VERBOSITY;
extern Bool BASE_MODULE_LOADED;

inline void checkConsole(SEXP ptr)
{
  // FIXME
}

void load_base_module();

template<typename StorageOrderType>
std::map<String, MultiArray> writeDataTable(SEXP data);


template<typename StorageOrderType>
SEXP readDataTable(const std::map<String, MultiArray> & dataMap);


IndexRange makeRange(const Rcpp::RObject & lower,
                            const Rcpp::RObject & upper);


template<typename StorageOrderType>
SEXP getMonitors(const std::map<String, NodeArrayMonitor> & monitorsMap, const String & type);

template <class InType> 
InType evalRfun(const Rcpp::Function & fun, const std::vector<InType> & invec) {

  InType outvec;
  switch(invec.size()) {
    case 1:
      outvec = fun(invec[0]);
      break;
    case 2:
      outvec = fun(invec[0], invec[1]);
      break;
    case 3:
      outvec = fun(invec[0], invec[1], invec[2]);
      break;
    case 4:
      outvec = fun(invec[0], invec[1], invec[2], invec[3]);
      break;
    case 5:
      outvec = fun(invec[0], invec[1], invec[2], invec[3], invec[4]);
      break;
    default:
      throw Biips::RuntimeError("Too many arguments in R function. must be <= 5");
      break;
      }
    return outvec;
}

Rcpp::NumericVector arrayToVector(const Biips::NumArray & array );


#endif /* RBIIPS_UTILS_H_ */
