library(stringr)


makelinegraph <- function(gtab) {
  # take a matrix of edges, return a matrix of edges (of the line graph)
  #
  linegraph <- data.frame()
  for (i in 1:nrow(gtab)) {
    for (j in 1:nrow(gtab)) {
      if (gtab[i,2] == gtab[j,1]) {
        a <- str_join(gtab[i,1],gtab[i,2],sep="-")
        b <- str_join(gtab[j,1],gtab[j,2],sep="-")
        c <- gtab[i,3] + gtab[j,3]
        linegraph <- rbind(linegraph, data.frame(From=a, To=b, Wt=c))
      }
    }
  }
  linegraph
}
