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
#' 
#' 
#' 
#' 
#' compareMtoK(handy, hmat(c(2,2), as.list(1:2)))
#' compareMtoK(handy, hmat(c(3,3), as.list(1:2)))
#' compareMtoK(handy, hmat(c(4,4,4), as.list(1:3)))
#' compareMtoK(handy, hmat(c(3,3,2,2), as.list(1:4)))
library(algstat)
options(max.print=1000000)

doMetropolis <- function(init, mat) {
  
  out <- metropolis(tab2vec(init), mat, 1e6, engine = "Cpp")
  str(out)
}

compareMtoK <- function(init, mat) {
  #moves <- markov(mat)
  cat("\nMARKOV-------------------\n")
  #print(moves)
  #print(dim(moves))
  #doMetropolis(init, moves)
  cat("\n\nKERNEL-------------------\n")
  kernMoves <- getKernMoves(mat)
#  print(kernMoves)
  print(dim(kernMoves))
  #doMetropolis(init, moves)
}
