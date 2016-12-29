#include "rbiips_utils.h"
#include <compiler/Compiler.hpp>
#include <BiipsBase.hpp>


Size VERBOSITY = 1;
Bool BASE_MODULE_LOADED = false;


void load_base_module()
{
  using namespace Biips;
  BEGIN_RBIIPS
  loadBaseModule(Compiler::FuncTab(), Compiler::DistTab());
  BASE_MODULE_LOADED = true;
  VOID_END_RBIIPS
}


template<>
std::map<String, MultiArray> writeDataTable<ColumnMajorOrder>(SEXP data)
{
  std::map<String, MultiArray> data_map;

  if (VERBOSITY>1)
    rbiips_cout << PROMPT_STRING << "Writing data table" << endl;

  Rcpp::List data_list(data);
  if (!data_list.hasAttribute("names"))
  {
    rbiips_cerr << "Warning: Missing variable names" << endl;
    return data_map;
  }

  if (VERBOSITY>1)
    rbiips_cout << INDENT_STRING << "Variables:";

  Rcpp::CharacterVector names = data_list.attr("names");
  for (int i=0; i<names.size(); ++i)
  {
    String var_name(names[i]);
    if (VERBOSITY>1)
      rbiips_cout << " " << var_name;

    Rcpp::NumericVector r_vec = data_list[var_name];
    MultiArray marray;

    if (!r_vec.hasAttribute("dim"))
    {
      DimArray::Ptr p_dim(new DimArray(1, r_vec.size()));
      ValArray::Ptr p_val(new ValArray(r_vec.size()));
      std::replace_copy(r_vec.begin(), r_vec.end(), p_val->begin(), NA_REAL, BIIPS_REALNA);
      marray.SetPtr(p_dim, p_val);
    }
    else
    {
      Rcpp::IntegerVector r_dim = r_vec.attr("dim");
      DimArray::Ptr p_dim(new DimArray(r_dim.begin(), r_dim.end()));
      ValArray::Ptr p_val(new ValArray(r_vec.size()));
      std::replace_copy(r_vec.begin(), r_vec.end(), p_val->begin(), NA_REAL, BIIPS_REALNA);
      marray.SetPtr(p_dim, p_val);
    }

    data_map[var_name] = marray;
  }
  if (VERBOSITY>1)
    rbiips_cout << endl;

  return data_map;
}


template<>
SEXP readDataTable<ColumnMajorOrder>(const std::map<String, MultiArray> & dataMap)
{
  if (VERBOSITY>1)
    rbiips_cout << PROMPT_STRING << "Reading data table" << endl;

  Rcpp::List data_list;

  if (VERBOSITY>1)
    rbiips_cout << INDENT_STRING << "Variables:";

  Rcpp::CharacterVector names;
  std::map<String, MultiArray>::const_iterator it_table = dataMap.begin();
  for (; it_table!=dataMap.end(); ++it_table)
  {
    const String & var_name = it_table->first;
    const MultiArray & values_array = it_table->second;

    // dim
    Rcpp::IntegerVector dim(values_array.Dim().begin(), values_array.Dim().end());

    Size len = values_array.Dim().Length();
    Rcpp::NumericVector values(len);

    std::replace_copy(values_array.Values().begin(), values_array.Values().end(), values.begin(), BIIPS_REALNA, NA_REAL);

    values.attr("dim") = dim;

    data_list[var_name] = values;

    if (VERBOSITY>1)
      rbiips_cout << " " << var_name;
  }
  if (VERBOSITY>1)
    rbiips_cout << endl;

  return data_list;
}


IndexRange makeRange(const Rcpp::RObject & lower,
                            const Rcpp::RObject & upper)
{
  if (lower.isNULL() || upper.isNULL())
    return IndexRange();

  Rcpp::IntegerVector il(lower);
  Rcpp::IntegerVector iu(upper);
  if (il.size() != iu.size())
    throw LogicError("length mismatch between lower and upper limits");

  IndexRange::Indices lind(il.begin(), il.end());
  IndexRange::Indices uind(iu.begin(), iu.end());

  IndexRange r = IndexRange(lind, uind);
  return r;
}


template<>
SEXP getMonitors<ColumnMajorOrder>(const std::map<String, NodeArrayMonitor> & monitorsMap, const String & type)
{
  Rcpp::List smcarray_list;

  std::map<String, NodeArrayMonitor>::const_iterator it_map;
  for (it_map = monitorsMap.begin(); it_map != monitorsMap.end(); ++it_map)
  {
    const String & name = it_map->first;
    const NodeArrayMonitor & monitor = it_map->second;

    // dim
    Rcpp::IntegerVector dim_particles(monitor.GetValues().Dim().begin(), monitor.GetValues().Dim().end());
    Rcpp::IntegerVector dim_array(monitor.GetRange().Dim().begin(), monitor.GetRange().Dim().end());

    // names(dim)
    Rcpp::CharacterVector dim_names(dim_particles.size(), "");
    dim_names[dim_names.size()-1] = "particle";

    dim_particles.attr("names") = dim_names;

    Size len_part = monitor.GetValues().Dim().Length();
    Rcpp::NumericVector values(len_part);
    const ValArray & values_val = monitor.GetValues().Values();
    std::replace_copy(values_val.begin(), values_val.end(), values.begin(), BIIPS_REALNA, NA_REAL);
    values.attr("dim") = dim_particles;

    const ValArray & weight_val = monitor.GetWeights().Values();
    Rcpp::NumericVector weights(weight_val.begin(), weight_val.end());
    weights.attr("dim") = dim_particles;

    const ValArray & ess_val(monitor.GetESS().Values());
    Rcpp::NumericVector ess(ess_val.begin(), ess_val.end());
    ess.attr("dim") = dim_array;

    const ValArray & discrete_val(monitor.GetDiscrete().Values());
    Rcpp::LogicalVector discrete(discrete_val.begin(), discrete_val.end());
    discrete.attr("dim") = dim_array;

    const ValArray & iter_val(monitor.GetIterations().Values()+1);
    Rcpp::NumericVector iterations(iter_val.begin(), iter_val.end());
    iterations.attr("dim") = dim_array;

    const Types<Types<String>::Array>::Array & cond = monitor.GetConditionalNodeNames();
    Size len = monitor.GetRange().Length();
    Rcpp::List cond_list(len);
    Rcpp::CharacterVector cond_vec;
    if (cond.size() == len) {
      for (Size i=0; i < len; ++i)
      {
        cond_list[i] = Rcpp::CharacterVector(cond[i].begin(), cond[i].end());
      }
      cond_list.attr("dim") = dim_array;
    }
    else if (cond.size() == 1) {
      cond_vec.assign(cond[0].begin(), cond[0].end());
    }
    else {
      throw LogicError("conditionals must either be of the same size as the node array or of size 1.");
    }


    const IndexRange::Indices & lower_ind = monitor.GetRange().Lower();
    Rcpp::IntegerVector lower(lower_ind.begin(), lower_ind.end());

    const IndexRange::Indices & upper_ind = monitor.GetRange().Upper();
    Rcpp::IntegerVector upper(upper_ind.begin(), upper_ind.end());

    Rcpp::List smcarray;
    smcarray["values"] = values;
    smcarray["weights"] = weights;
    smcarray["ess"] = ess;
    smcarray["discrete"] = discrete;
    smcarray["iterations"] = iterations;
    if (cond.size() == len)
      smcarray["conditionals"] = cond_list;
    else
      smcarray["conditionals"] = cond_vec;
    smcarray["name"] = Rcpp::wrap(monitor.GetName());
    smcarray["lower"] = lower;
    smcarray["upper"] = upper;
    smcarray["type"] = Rcpp::wrap(type);

    smcarray.attr("class") = "smcarray";

    smcarray_list[name] = smcarray;
  }

  return smcarray_list;
}


Rcpp::NumericVector arrayToVector(const Biips::NumArray & array ) {
  Rcpp::NumericVector vec(array.Values().begin(), array.Values().end());
  Rcpp::IntegerVector dim(array.Dim().begin(), array.Dim().end());
  vec.attr("dim") = dim;
  return vec;
}
