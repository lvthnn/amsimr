#ifndef AMSIMR_BINDINGS_UTILS_H
#define AMSIMR_BINDINGS_UTILS_H

// [[Rcpp::plugins(cpp20))]]
#include <Rcpp.h>

#include <amsim/component_type.h>
#include <amsim/mating_type.h>
#include <amsim/log_level.h>

#include <optional>
#include <string>

template<typename T>
inline std::optional<T> _optional(SEXP x) {
  if (x == R_NilValue) return std::nullopt;
  return Rcpp::as<T>(x);
}

inline amsim::MatingType _s_MatingType(const std::string& mating_type) {
  if (mating_type == "random")
    return amsim::MatingType::RANDOM;
  else if (mating_type == "assortative")
    return amsim::MatingType::ASSORTATIVE;
  else
    Rcpp::stop("Unknown MatingType " + mating_type);
}

inline amsim::ComponentType _s_ComponentType(const std::string& component_type) {
  if (component_type == "genetic")
    return amsim::ComponentType::GENETIC;
  else if (component_type == "environmental")
    return amsim::ComponentType::ENVIRONMENTAL;
  else if (component_type == "vertical")
    return amsim::ComponentType::VERTICAL;
  else if (component_type == "total")
    return amsim::ComponentType::TOTAL;
  else
    Rcpp::stop("Unknown ComponentType " + component_type);
}

inline amsim::LogLevel _s_LogLevel(const std::string& log_level) {
  return amsim::LogLevel::DEBUG;
}

#endif // AMSIMR_BINDINGS_UTILS_H
