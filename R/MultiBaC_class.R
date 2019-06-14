#' Class MultiBaC
#' @return Instance of class "plsClass"
#' #' \describe{
#'     \item{X}{A list containing pls results of matrix X}
#'     \item{Y}{A list containing pls results of matrix Y}
#'     \item{Bhat}{Coefficients matrix of the PLS model}
#'     \item{scale}{TRUE or FALSE. Whether X and Y have been scaled or not}
#'     \item{center}{TRUE or FALSE. Whether X and Y have been centered or not}
#' }
multibacClass = function(multibac, ...) {
  model = structure(multibac, class = "multibacClass",
                    my_attribute = "MultiBaC object")
  return(model)
}


