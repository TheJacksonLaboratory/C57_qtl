read_neogen <- function(inputfiles, outfile, B6J=c("C57BL/6J_A", "C57BL/6J_B"), B6N=c("C57L/J_A", "C57L/J_B")){
  # Read Karl Broman analysis result with unique and mapped markers, save a list of the unique only
  kb <- read.csv(url("https://raw.githubusercontent.com/kbroman/MUGAarrays/master/UWisc/mini_uwisc_v2.csv"))
  kblist <- kb$marker[kb$unique & (!kb$unmapped)]
  # Read the Neogen input file and write the genotypes matrix into the output csv file
  tbl <- NULL
  for (inputfile in inputfiles){
    tbl1 <- read.delim(inputfile, skip = 10, header = F)
    tbl <- rbind(tbl, tbl1)
  }
  w3 <- reshape(tbl[,1:3], v.names="V3", idvar = "V1", timevar="V2", direction="wide")
  rownames(w3) <- w3$V1
  w3 <- w3[,-1]
  colnames(w3) <- gsub("^V3.", "", colnames(w3))
  w3 <- t(w3)
  w4 <- reshape(tbl[,c(1:2, 4)], v.names="V4", idvar = "V1", timevar="V2", direction="wide")
  rownames(w4) <- w4$V1
  w4 <- w4[,-1]
  colnames(w4) <- gsub("^V4.", "", colnames(w4))
  w4 <- t(w4)
  
  # Find the B6J and B6N consensus. If one allele is different remove the marker
  # If B6J==B6N remove the marker
  keepr <- apply(w3[B6J,], 2, function(col) length(unique(col)) == 1) & 
    apply(w4[B6J,], 2, function(col) length(unique(col)) == 1) &
    apply(w3[B6N,], 2, function(col) length(unique(col)) == 1) &
    apply(w4[B6N,], 2, function(col) length(unique(col)) == 1) &
    w3[B6J[1],] == w4[B6J[1],] &
    w3[B6N[1],] == w4[B6N[1],] &
    w3[B6J[1],] != w3[B6N[1],]
  
  # Set the allele names to A and B for B6J and B6N
  w3 <- w3[, keepr]
  w4 <- w4[, keepr]
  d3 <- w3 == rep(w3[B6J[1],],each=nrow(w3))
  d3[d3==TRUE] = "A"
  d3[d3==FALSE] = "B"
  d4 <- w4 == rep(w4[B6J[1],],each=nrow(w4))
  d4[d4==TRUE] = "A"
  d4[d4==FALSE] = "B"
  
  # Concatenate the two alleles
  for (i in 1:nrow(d3)) d3[i,] <- paste0(d3[i,], d4[i,])
  d3 <- d3[, intersect(colnames(d3), kblist)]
  # Separate the founders and write them in the first two rows followed by the rest
  dout <- rbind(d3[B6J[1],,drop=F], d3[B6N[1],,drop=F], d3[!rownames(d3) %in% c(B6N, B6J),])
  rownames(dout) <- c("C57BL_6J", "C57L_6J", rownames(dout)[3:nrow(dout)])
  write.csv(dout, file=outfile)
} 
