#' @examples
#' # From algstat documentation
#' data(handy)
#' 
#' exp   <- loglin(handy, as.list(1:2), fit = TRUE)$fit
#' e <- unname(tab2vec(exp))
#' h <- t(t(unname(tab2vec(handy))))
#' chisq <- algstat:::computeX2sCpp(h, e)
#' 
#' out <- loglinear(~ Gender + Handedness, data = handy)
#' chisqs <- algstat:::computeX2sCpp(out$steps, e)
#' 
#' mean(chisqs >= chisq)
#' fisher.test(handy)$p.value
#' 
#' 
#' 
#' 
#' 
#' A <- hmat(c(2,2), as.list(1:2))
#' moves <- markov(A)
#' outMarkov <- metropolis(tab2vec(handy), moves, 1e4, engine = "Cpp")
#' str(outMarkov)
#' 
#' kernMoves <- matrix(m2_kernel(m2_matrix(A)))
#' outMarkov <- metropolis(tab2vec(handy), kernMoves, 1e4, engine = "Cpp")
#' 
#' # showSteps(out$steps)
#' 
#' 
#' library(microbenchmark)
#' microbenchmark(
#'   metropolis(tab2vec(handy), moves, engine = "Cpp"),
#'   metropolis(tab2vec(handy), moves, engine = "R")
#' )
#' 
#' # cpp ~ 20-25x faster
library(m2r)
library(algstat)

doMetropolis <- function(init, mat) {
  
  out <- metropolis(tab2vec(init), mat, 1e6, engine = "Cpp")
  str(out)
}

compareMtoK <- function(init, mat) {
  moves <- markov(mat)
  cat("\nMARKOV-------------------")
  doMetropolis(init, moves)
  cat("\n\nKERNEL-------------------")
  kernMoves <- matrix(m2_kernel(m2_matrix(A)))
  doMetropolis(init, moves)
}
compareMtoK(handy, A <- hmat(c(2,2), as.list(1:2)))
