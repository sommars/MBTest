make_matrix_str <- function(mat) {
  matrix_str <- "matrix{"
  for (i in 1:nrow(mat)) {
    matrix_str <- paste0(matrix_str,"{",paste0(mat[i,],collapse=","), "}")
    if (i != nrow(mat)) {
      matrix_str <- paste0(matrix_str, ",")
    } else {
      matrix_str <- paste0(matrix_str, "}")
    }
  }
  paste("kernel", matrix_str)
}

parse_kernel <- function(kern) {
  # Could be better
  splitKern <- strsplit(kern, "\\^")
  for (i in 1:nchar(splitKern[[1]][2])) {
    if (substr(splitKern[[1]][2],i,i) == ",") {
      rowCount <- as.numeric(substr(splitKern[[1]][2],1,i-1))
      break
    }
  }
  
  for (i in 1:nchar(splitKern[[1]][3])) {
    if (substr(splitKern[[1]][3],i,i) == ",") {
      colCount <- as.numeric(substr(splitKern[[1]][3],1,i-1))
      matrix_str <- substr(splitKern[[1]][3],i+1,nchar(splitKern[[1]][3])-1)
      matrix_str <- str_replace_all(str_replace_all(matrix_str, "\\{", ""),"\\}","")
      break
    }
  }
  matrix(lapply(strsplit(matrix_str,","), as.numeric)[[1]], nrow = rowCount, ncol = colCount, byrow = TRUE)
}

run_m2 <- function(code) {
  dir <- getwd()
  codeFile <- "M2FileCreatedByR.M2"
  outFile <- "M2OutputForR"
  
  code <- paste('"', outFile, '"<< toExternalString(', code, ") << endl << close", sep = "")
  writeLines(code, con = paste(dir, codeFile, sep = "/"))
  system2("M2", paste("--script", codeFile))
  m2_output <- readLines(paste(dir, outFile, sep = "/"))
  system2("rm", codeFile)
  system2("rm", outFile)
  m2_output
}

getKernMoves <- function(mat) {
  parse_kernel(run_m2(make_matrix_str(mat)))
}
