summary.KScens <- function(object, outp = c("list", "table"),
                            print.AIC = TRUE, print.BIC = TRUE,
                            print.infoBoot = FALSE, ...) {
  if (!inherits(object, "KScens")) {
    stop("Use only 'KScens' objects")
  }
  outp <- match.arg(outp)
  if (!outp %in% c("list", "table")) {
    stop("Invalid value of outp. Use 'table' or 'list'.")
  }
  object$outp <- outp
  object$print.AIC <- print.AIC
  object$print.BIC <- print.BIC
  object$print.infoBoot <- print.infoBoot
  class(object) <- c("summary.KScens", class(object))
  object
}
