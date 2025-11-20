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
