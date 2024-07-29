print.CvMcens <- function(x, prnt = TRUE, outp = c("list", "table"), ...) {
  outp <- match.arg(outp)
  if (prnt && !outp %in% c("list", "table")) {
    stop("Invalid value of outp. Use 'table' or 'list'.")
  }
  if (prnt) {
    if (outp == "table") {
      cat("Distribution:", x$Distribution, "\n")
      if (!is.null(x$Hypothesis)) {
        cat("\nNull hypothesis:\n")
        header1 <- c("Parameter", "Value")
        max_col_width1 <- max(nchar(header1), nchar(names(x$Hypothesis)))
        cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                    strrep("-", max_col_width1)))
        cat(sprintf("%-*s | %-*s\n", max_col_width1, header1[1],
                    max_col_width1, header1[2]))
        cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                    strrep("-", max_col_width1)))
        for (i in 1:length(x$Hypothesis)) {
          cat(sprintf("%-*s | %-*s\n", max_col_width1, names(x$Hypothesis)[i],
                      max_col_width1, unname(x$Hypothesis)[i]))
        }
        cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                    strrep("-", max_col_width1)))
      }
      cat("\nCvM Test results:\n")
      header <- c("Metric", "Value")
      max_col_width <- max(nchar(header), nchar(names(x$Test)))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                  strrep("-", max_col_width)))
      cat(sprintf("%-*s | %-*s\n", max_col_width, header[1],
                  max_col_width, header[2]))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                  strrep("-", max_col_width)))
      for (i in 1:length(x$Test)) {
        cat(sprintf("%-*s | %-*s\n", max_col_width, names(x$Test)[i],
                    max_col_width, unname(x$Test)[i]))
      }
      cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                  strrep("-", max_col_width)))
      cat("\nParameter estimates:\n")
      header1 <- c("Parameter", "Value")
      max_col_width1 <- max(nchar(header1), nchar(names(x$Estimates)))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      cat(sprintf("%-*s | %-*s\n", max_col_width1, header1[1],
                  max_col_width1, header1[2]))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      for (i in 1:length(x$Estimates)) {
        cat(sprintf("%-*s | %-*s\n", max_col_width1, names(x$Estimates)[i],
                    max_col_width1, unname(x$Estimates)[i]))
      }
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
    } else {
      cat("Distribution:", x$Distribution, "\n")
      if (!is.null(x$Hypothesis)) {
        cat("\nNull hypothesis:\n")
        print(x$Hypothesis)
      }
      cat("\nCvM Test results:\n")
      print(x$Test)
      cat("\nParameter estimates:\n")
      print(x$Estimates)
    }
  }
  invisible(x)
}
