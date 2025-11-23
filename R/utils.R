#' Assert valid probability values
#'
#' @description
#' Validates that a numeric vector contains valid probability values between
#' 0 and 1 with no missing values.
#'
#' @param x Numeric vector to validate.
#' @param len Integer or NULL. Expected length of the vector (default: NULL,
#'   no length check).
#'
#' @return Invisibly returns `x` if validation passes; otherwise throws an
#'   error.
#'
#' @keywords internal
assert_probs <- function(x, len = NULL) {
  checkmate::assert_numeric(
    x,
    len = len,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    all.missing = FALSE
  )
}

#' Assert valid correlation values
#'
#' @description
#' Validates that a numeric vector contains valid correlation values between
#' -1 and 1 with no missing values.
#'
#' @param x Numeric vector to validate.
#' @param len Integer or NULL. Expected length of the vector (default: NULL,
#'   no length check).
#'
#' @return Invisibly returns `x` if validation passes; otherwise throws an
#'   error.
#'
#' @keywords internal
assert_cors <- function(x, len = NULL) {
  checkmate::assert_numeric(
    x,
    len = len,
    lower = -1,
    upper = 1,
    any.missing = FALSE,
    all.missing = FALSE
  )
}
