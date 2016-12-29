#include "rbiips_utils.h"
#include <common/Accumulator.hpp>
#include <common/Utility.hpp>

using namespace Biips;
using std::endl;

using namespace Biips;

template <typename FeatureIterator>
static Accumulator accumulate(SEXP values, SEXP weights, FeatureIterator firstFeat, FeatureIterator lastFeat)
{
  Rcpp::NumericVector values_vec(values);
  Rcpp::NumericVector weights_vec(weights);
  if (values_vec.size() != weights_vec.size())
    throw LogicError("values and weights must have same length.");

  Accumulator accu;
  for (;firstFeat != lastFeat; ++firstFeat)
    accu.AddFeature(*firstFeat);

  accu.Init();

  for (int i = 0; i<values_vec.size(); ++i)
    accu.Push(values_vec[i], weights_vec[i]);

  return accu;
}


RcppExport SEXP wtd_stat (SEXP values, SEXP weights, SEXP order)
{
  BEGIN_RBIIPS

  Size order_s = Rcpp::as<Size>(order);

  static StatTag features[] = {MEAN, VARIANCE, SKEWNESS, KURTOSIS};
  Accumulator accu = accumulate(values, weights, features, features + order_s);

  Rcpp::NumericVector stat_vec(order_s);

  stat_vec[0] = accu.Mean();
  if (order_s>=2)
    stat_vec[1] = accu.Variance();
  if (order_s>=3)
    stat_vec[2] = accu.Skewness();
  if (order_s>=4)
    stat_vec[3] = accu.Kurtosis();

  static String names[] = {"mean", "var", "skew", "kurt"};

  Rcpp::CharacterVector names_vec(names, names + order_s);

  stat_vec.attr("names") = names_vec;

  return stat_vec;
  END_RBIIPS
}


RcppExport SEXP wtd_mean (SEXP values, SEXP weights)
{
  BEGIN_RBIIPS

  static StatTag features[] = {MEAN};
  Accumulator accu = accumulate(values, weights, features, features + sizeof(features)/sizeof(StatTag));

  return Rcpp::wrap(accu.Mean());
  END_RBIIPS
}


RcppExport SEXP wtd_var (SEXP values, SEXP weights)
{
  BEGIN_RBIIPS

  static StatTag features[] = {VARIANCE};
  Accumulator accu = accumulate(values, weights, features, features + sizeof(features)/sizeof(StatTag));

  return Rcpp::wrap(accu.Variance());
  END_RBIIPS
}


RcppExport SEXP wtd_skew (SEXP values, SEXP weights)
{
  BEGIN_RBIIPS

  static StatTag features[] = {SKEWNESS};
  Accumulator accu = accumulate(values, weights, features, features + sizeof(features)/sizeof(StatTag));

  return Rcpp::wrap(accu.Skewness());
  END_RBIIPS
}


RcppExport SEXP wtd_kurt (SEXP values, SEXP weights)
{
  BEGIN_RBIIPS

  static StatTag features[] = {KURTOSIS};
  Accumulator accu = accumulate(values, weights, features, features + sizeof(features)/sizeof(StatTag));

  return Rcpp::wrap(accu.Kurtosis());
  END_RBIIPS
}


RcppExport SEXP wtd_quantile (SEXP values, SEXP weights, SEXP probs)
{
  BEGIN_RBIIPS

  Rcpp::NumericVector values_vec(values);
  Rcpp::NumericVector weights_vec(weights);
  if (values_vec.size() != weights_vec.size())
    throw LogicError("values and weights must have same length.");
  Rcpp::NumericVector probs_vec(probs);
  Rcpp::NumericVector quant_vec(probs_vec.size());

  QuantileAccumulator accu(probs_vec.begin(), probs_vec.end());
  accu.Init();

  for (int i = 0; i<values_vec.size(); ++i)
    accu.Push(values_vec[i], weights_vec[i]);

  for (int i = 0; i<probs_vec.size(); ++i)
    quant_vec[i] = accu.Quantile(i);

  quant_vec.attr("names") = Rcpp::CharacterVector(probs);

  return quant_vec;
  END_RBIIPS
}


RcppExport SEXP wtd_median (SEXP values, SEXP weights)
{
  BEGIN_RBIIPS

  Rcpp::NumericVector values_vec(values);
  Rcpp::NumericVector weights_vec(weights);
  if (values_vec.size() != weights_vec.size())
    throw LogicError("values and weights must have same length.");

  std::vector<double> proba(1, 0.5);

  QuantileAccumulator accu(proba.begin(), proba.end());
  accu.Init();

  for (int i = 0; i<values_vec.size(); ++i)
    accu.Push(values_vec[i], weights_vec[i]);

  return Rcpp::wrap(accu.Quantile(0));
  END_RBIIPS
}



RcppExport SEXP wtd_table(SEXP values, SEXP weights)
{
  BEGIN_RBIIPS

  Rcpp::NumericVector values_vec(values);
  Rcpp::NumericVector weights_vec(weights);
  if (values_vec.size() != weights_vec.size())
    throw LogicError("values and weights must have same length.");

  DiscreteAccumulator accu;
  accu.Init();

  for (int i = 0; i<values_vec.size(); ++i)
    accu.Push(values_vec[i], weights_vec[i]);

  const DiscreteHistogram & hist = accu.Pdf();
  Types<Scalar>::Array vec = hist.GetPositions();
  Rcpp::CharacterVector x(vec.begin(), vec.end());
  vec = hist.GetFrequencies();
  Rcpp::NumericVector table(vec.begin(), vec.end());
  table.attr("names") = x;
  table.attr("class") = "table";

  return table;
  END_RBIIPS
}


RcppExport SEXP wtd_mode(SEXP values, SEXP weights)
{
  BEGIN_RBIIPS

  Rcpp::NumericVector values_vec(values);
  Rcpp::NumericVector weights_vec(weights);
  if (values_vec.size() != weights_vec.size())
    throw LogicError("values and weights must have same length.");

  DiscreteAccumulator accu;
  accu.Init();

  for (int i = 0; i<values_vec.size(); ++i)
    accu.Push(values_vec[i], weights_vec[i]);

  return Rcpp::wrap(accu.Mode());
  END_RBIIPS
}
