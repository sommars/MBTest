#-------------------------------------------------------------------------------
do_scrape <- function(df) {
  dir <- getwd()
  outFile <- paste(dir, "ScrapeOutput.txt", sep = "/")
  hmatFileLocation <- paste(dir, "RHMatFile", sep = "/")
  mbFileLocation <- paste(dir, "RMBFile", sep = "/")
  
  for (dfRow in 1:dim(df)[1]) {
    testRow <- df[dfRow,]
    fileName <- testRow$name
    print(fileName)
    write("NEWSYSTEM",file = outFile,append = TRUE)
    write(as.character(fileName),file = outFile,append = TRUE)
    write(paste(c("Degree: ", testRow$deg),collapse = ""),file = outFile,append = TRUE)
    
    hmatFile <- sprintf("http://markov-bases.de/data/%s/%s.mat", fileName, fileName)
    basisFile <- sprintf("http://markov-bases.de/data/%s/%s.mar", fileName, fileName)
    
    if (download.file(hmatFile, destfile = hmatFileLocation, quiet = TRUE) != 0)
      next
    myHmat <- as.matrix(read.table(hmatFileLocation, sep = " ", as.is=TRUE, skip = 1, strip.white = TRUE))
    myHmat <- myHmat[,-dim(myHmat)[[2]]]
    kernMoves <- getKernMoves(myHmat)
    if (download.file(basisFile, destfile = mbFileLocation, quiet = TRUE) != 0)
      next
    cleanedUpFile <- gsub( "  ", " ", readLines(mbFileLocation))
    cat(cleanedUpFile, file=mbFileLocation, sep="\n")
    markMoves <- as.matrix(read.table(mbFileLocation, sep = c(" ", "  "), as.is=TRUE, skip = 1, strip.white = TRUE))
    markNrows <- dim(markMoves)[[1]]
    markMoves <- matrix(markMoves[,-dim(markMoves)[[2]]], nrow = markNrows)
    markMoves <- matrix(markMoves[,-1], nrow = markNrows)
    markMoves <- t(markMoves)
    
    markMovesRows <- dim(markMoves)[[1]]
    markMoves <- matrix(markMoves[, !duplicated(t(markMoves))], nrow = markMovesRows)
    
    kernMoves <- cbind(kernMoves, kernMoves * -1)
    f.obj <- vector(mode="numeric", length = dim(kernMoves)[[2]])
    f.con <- kernMoves
    f.dir <- vector(mode="character", length = dim(kernMoves)[[1]])
    for (i in 1:length(f.dir)) {
      f.dir[i] <- "="
    }
    lpResults <- c()
    for (i in 1:dim(markMoves)[[2]]) {
      lpResult <- lp ("max", f.obj, f.con, f.dir, c(markMoves[,i]), all.int = TRUE)
      write(paste(lpResult$solution,collapse = " "),file = outFile,append = TRUE)
      if (sum(lpResult$solution) > testRow$deg) {
        write("FirstGuessWrong",file = outFile,append = TRUE)
      }
      lpResults <- c(lpResults, sum(lpResult$solution))
    }
    write(paste(c("Max value: ", max(lpResults)),collapse = ""),file = outFile,append = TRUE)
    write(paste(c("Min value: ", min(lpResults)),collapse = ""),file = outFile,append = TRUE)
    write(paste(c("Median value: ", median(lpResults)),collapse = ""),file = outFile,append = TRUE)
  }
}

#-------------------------------------------------------------------------------
read_data_file <- function() {
  dir <- getwd()
  file <- "DBDump.txt"
  read.csv2(paste(dir, file, sep = "/"), sep = '"')
}


#df <- read_data_file()
#do_scrape(df[9,])
