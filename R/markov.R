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
#' help(loglin) ... HairEyeColor
#' data(handy)
#' compareMtoK(handy, hmat(c(4,4), as.list(1:3)))
#' hec <- HairEyeColor
#' compareMtoK(HairEyeColor, hmat(c(3,3,2), list(c(1,2),c(1,3),c(2,3))))
#' compareMtoK(HairEyeColor, hmat(c(4,4,2), list(c(1,2),c(1,3),c(2,3))))
#' compareMtoK(handy, hmat(c(7,4,4), as.list(1:3)))
#' compareMtoK(handy, hmat(c(5,10,6,5), as.list(1:4)))
#'
#' ---------------------------
#' 3e6 for c(3,3,2)
#' Unique in Markov:  1477 
#' Unique in kernel:  1466 
#' In both:  8688 
#'
#'
#' set_latte_path("/home/jeff/Desktop/Software/latte-integrale-1.7.3b/dest/bin")
#' data(politics)
#' countTables(politics)
#' (A <- hmat(c(2,2), list(1, 2)))
#' countTables(politics, A)
#'
#' data(HairEyeColor)
#' df <- as.data.frame(HairEyeColor)
#' aa <- df[df$Hair != "Blond" & df$Eye != "Green",]
#' aa$Hair <- factor(aa$Hair)
#' aa$Eye <- factor(aa$Eye)
#' tab <- xtabs(Freq ~ ., aa)
#' compareMtoK(tab, hmat(c(3,3,2), list(c(1,2),c(1,3),c(2,3))))
#' compareMtoK(HairEyeColor, hmat(c(4,4,2), list(c(1,2),c(1,3),c(2,3))))
#' 
#' 
#' compareMtoK(HairEyeColor, hmat(c(4,3,2,3), list(c(1,2),c(1,3),c(2,3),c(1,4))))
#' 
#' init <- tab
#' mat <- hmat(c(3,3,2), list(c(1,2),c(1,3),c(2,3)))
#' 
#' # countTables(tab, hmat(c(3,3,2), list(c(1,2),c(1,3),c(2,3))))   87,024
#' # countTables(HairEyeColor, hmat(c(4,4,2), list(c(1,2),c(1,3),c(2,3))))  10,822,861,056
library(algstat)
library(lpSolve)
library(rvest)
library(XML)
data(handy)


options(max.print=1000000)

#-------------------------------------------------------------------------------
doMetropolis <- function(init, mat) {
  metropolis(tab2vec(init), mat, 1e6, engine = "Cpp")
}

#-------------------------------------------------------------------------------
compareMtoK <- function(init, mat) {
  cat("\nMARKOV-------------------\n")
  moves <- markov(mat)
  print(dim(moves))
  #markovMet <- metropolis(tab2vec(init), moves, 5e6, engine = "Cpp")
  cat("\n\nKERNEL-------------------\n")
  kernMoves <- getKernMoves(mat)
  print(dim(kernMoves))
  doCompare2(moves, kernMoves, init)
  #kernMet <- metropolis(tab2vec(init), kernMoves, 1e6, engine = "Cpp")
  #doCompare(markovMet$steps, kernMet$steps)
  #markSteps<-markovMet$steps
}

#-------------------------------------------------------------------------------
doCompare2 <- function(markSteps, kernMoves, init) {
  #markSteps <- markMoves
  init <- tab2vec(init)
  markStepRows <- dim(markSteps)[[1]]
  markSteps <- matrix(markSteps[, !duplicated(t(markSteps))], nrow = markStepRows)
  markOnly <- 0
  Both <- 0
  
  kernMoves <- cbind(kernMoves, kernMoves * -1)
  f.obj <- vector(mode="numeric", length = dim(kernMoves)[[2]])
  f.con <- kernMoves
  f.dir <- vector(mode="character", length = dim(kernMoves)[[1]])
  for (i in 1:length(f.dir)) {
    f.dir[i] <- "="
  }
  
  for (i in 1:dim(markSteps)[[2]]) {
    print(markSteps[,i])
    lpResult <- lp ("max", f.obj, f.con, f.dir, c(markSteps[,i]), all.int = TRUE)
    print(lpResult$solution)
    if (lpResult$status == 0) {
      Both <- Both + 1
    } else {
      markOnly <- markOnly + 1
    }
  }
  cat("Unique in Markov: ", markOnly, "\n")
  cat("In both: ", Both, "\n")
  if (markOnly != 0)
    stop("Darn")
  invisible(0)
}

#-------------------------------------------------------------------------------
cleanupMatrix <- function(M) {
  # Sort and remove duplicates
  M <- M[,do.call(order, as.data.frame(t(M))[1:dim(M)[[1]]])]
  M[, !duplicated(t(M))]
}

#-------------------------------------------------------------------------------
doCompare <- function(M1, M2) {
  #If M1 and M2 each find the same column twice, we'll record that 
  #only as 1. Not sure if this is best.
  M1Only <- 0
  M2Only <- 0
  Both <- 0
  M1 <- cleanupMatrix(M1)
  M2 <- cleanupMatrix(M2)
  M1Index <- 1
  M2Index <- 1
  RowCount <- dim(M1)[[1]]
  while (TRUE) {
    if (M1Index > dim(M1)[[2]]) {
      M2Only <- M2Only + (dim(M2)[[2]] - M2Index + 1)
      break
    }
    if (M2Index > dim(M2)[[2]]) {
      M1Only <- M1Only + (dim(M1)[[2]] - M1Index + 1)
      break
    }
    
    val <- "="
    for (i in 1:RowCount) {
      if (M1[,M1Index][i] < M2[,M2Index][i]) {
        val <- "M1<M2"
        break
      } else if (M1[,M1Index][i] > M2[,M2Index][i]) {
        val <- "M1>M2"
        break
      }        
    }
    if (val == "M1<M2") {
      M1Index <- M1Index + 1
      M1Only <- M1Only + 1
    } else if (val == "M1>M2") {
      M2Index <- M2Index + 1
      M2Only <- M2Only + 1
    } else if (val == "=") {
      currentCol <- M1[,M1Index]
      M1Index <- M1Index + 1
      M2Index <- M2Index + 1
      Both <- Both + 1
    }
  }  
  cat("Unique in Markov: ", M1Only, "\n")
  cat("Unique in kernel: ", M2Only, "\n")
  cat("In both: ", Both, "\n")
  invisible(0)
}
