#' Print method for class \code{"brar"}
#' @method print brar
#'
#' @param x Object of class \code{"brar"}
#' @param digits Number of digits for formatting of numbers
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"brar"} object
#'
#' @author Samuel Pawel
#'
#' @examples
#' brar_binomial(y1 = 20, y2 = 30, n1 = 50, n2 = 50, pH0 = 0.5)
#'
#' @noRd
#' @export
print.brar <- function(x, digits = 3, ...) {

    cat("DATA\n")
    print(x$dat, digits = digits, row.names = FALSE)
    cat("\n")

    cat("PRIOR PROBABILITIES\n")
    print(x$prior, digits = digits, row.names = FALSE)
    cat("\n")

    cat("BAYES FACTORS (BF_ij)\n")
    print(x$BF_ij, digits = digits, row.names = FALSE)
    cat("\n")

    cat("POSTERIOR PROBABILITIES\n")
    print(x$post, digits = digits, row.names = FALSE)
    cat("\n")

    cat("RANDOMIZATION PROBABILITIES\n")
    print(x$prand, digits = digits, row.names = FALSE)
    cat("\n")

    invisible(x)
}


## #' Plot method for class \code{"brar"}
## #' @method plot brar
## #'
## #' @param x Object of class \code{"brar"}
## #' @param ... Other arguments (for consistency with the generic)
## #'
## #' @return Plots
## #'
## #' @author Samuel Pawel
## #'
## #' @seealso \code{\link{brar}}
## #'
## #' @export
## plot.brar <- function(x, ...) {

## }
