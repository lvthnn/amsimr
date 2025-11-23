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
inline std::optional<T> optional(SEXP x) {
  if (x == R_NilValue) return std::nullopt;
  return Rcpp::as<T>(x);
}

inline amsim::MatingType strMatingType(const std::string& mating_type) {
  if (mating_type == "random")
    return amsim::MatingType::RANDOM;
  if (mating_type == "assortative")
    return amsim::MatingType::ASSORTATIVE;
  Rcpp::stop("Unknown MatingType " + mating_type);
}

inline amsim::ComponentType strComponentType(const std::string& component_type) {
  if (component_type == "genetic")
    return amsim::ComponentType::GENETIC;
  if (component_type == "environmental")
    return amsim::ComponentType::ENVIRONMENTAL;
  if (component_type == "vertical")
    return amsim::ComponentType::VERTICAL;
  if (component_type == "total")
    return amsim::ComponentType::TOTAL;
  Rcpp::stop("Unknown ComponentType " + component_type);
}

inline amsim::LogLevel strLogLevel(const std::string& log_level) {
  if (log_level == "debug")
    return amsim::LogLevel::DEBUG;
  if (log_level == "info")
    return amsim::LogLevel::INFO;
  if (log_level == "warning")
    return amsim::LogLevel::WARNING;
  if (log_level == "error")
    return amsim::LogLevel::ERROR;
  if (log_level == "NONE")
    return amsim::LogLevel::NONE;
  Rcpp::stop("Unknown LogLevel " + log_level);
}

#endif // AMSIMR_BINDINGS_UTILS_H
