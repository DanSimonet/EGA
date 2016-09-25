#'  Investigates the stability of EGA's estimation via bootstrap.
#'
#' \code{bootEGA} Estimates the number of dimensions of n bootstraps from the empirical correlation matrix,
#'  and returns a typical network (i.e. the network formed by the median pairwise partial correlations over the n bootstraps) and its dimensionality.
#'
#' @param data A dataframe with the variables to be used in the analysis
#' @param n An integer value representing the number of bootstraps
#'
#' @param medianStructure Logical. If true, returns the typical network of partial correlations (estimated via graphical lasso), which is the median of all pairwise correlations over the n bootstraps, and estimates its dimensions.
#' @param plot.MedianStructure Logical. If true, returns a plot of the typical network (partial correlations), which is the median of all pairwise correlations over the n bootstraps, and its estimated dimensions.
#' @author Hudson F. Golino <hfgolino at gmail.com>
#' @examples
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 500, medianStructure = TRUE, plot.MedianStructure = TRUE)
#' boot.intwl <- bootEGA(data = intelligenceBattery[,8:66], n = 500, medianStructure = TRUE, plot.MedianStructure = TRUE)
#'
#' \dontrun{
#' bootEGA(a)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export


#Bootstrap EGA:
bootEGA <- function(data, n, medianStructure = TRUE, plot.MedianStructure = TRUE){
  require(qgraph)
  require(bootnet)
  require(igraph)
  boot.ega <- bootnet(data, nBoot = n, default = "EBICglasso", computeCentrality = FALSE)
  bootGraphs <- vector("list", n)
  for (i in 1:n) {
    bootGraphs[[i]] <- boot.ega$boots[[i]]
  }
  boot.igraph <- vector("list", n)
  for (l in 1:n){
    boot.igraph[[l]] <-as.igraph(plot(bootGraphs[[l]], DoNotPlot = TRUE))
  }
  boot.wc <- vector("list", n)
  for (m in 1:n){
    boot.wc[[m]] <-walktrap.community(boot.igraph[[m]])
  }
  boot.ndim <- matrix(NA,nrow=i,ncol=2)
  for (m in 1:n){
    boot.ndim[m,2] <-max(boot.wc[[m]]$membership)
  }
  colnames(boot.ndim) <- c("Boot.Number", "N.Dim")
  boot.ndim[,1] <- seq_len(n)

  #Computing the median EBICglasso graph and its dimensionality:
  if(medianStructure == TRUE){
    mStructure <- list()
    for (k in 1:length(boot.ega$boots)) {
      mStructure[[k]]  <- boot.ega$boots[[k]][[1]]
      median.Structure <- apply(simplify2array(mStructure), 1:2, median)
    }
    median.igraph <- as.igraph(qgraph(median.Structure, DoNotPlot = TRUE))
    median.wc <-walktrap.community(median.igraph)
    median.ndim <- max(median.wc$membership)
    dim.variables <- data.frame(items = names(data), dimension = median.wc$membership)


  }

  #Ploting the mean EBICglasso graph:

  if (plot.MedianStructure == TRUE) {
    plot.median.ega <- qgraph(median.Structure, layout = "spring", vsize = 4, groups = as.factor(median.wc$membership))
  }


  #Computing the Median, SD, SE and CI of the estimated dimensions
  Median <- median(boot.ndim[,2])
  sd.boot <- sd(boot.ndim[,2])
  se.boot <- (1.253*sd.boot)/sqrt(nrow(boot.ndim))
  ciMult <- qt(.95/2 + .5, nrow(boot.ndim)-1)
  ci <- se.boot * ciMult
  summary.table <- data.frame(n.Boots = n, median.dim = Median, SD.dim = sd.boot, SE.dim = se.boot, CI.dim = ci, Lower = Median-ci, Upper = Median + ci)

  #Print results
  result <- list()
  result$n <- n
  result$boot.ndim <- boot.ndim
  result$bootGraphs <- bootGraphs
  result$summary.table <- summary.table
  medianGraph <- list()
  medianGraph$graph <- median.Structure
  medianGraph$median.dim.variables <- dim.variables[order(dim.variables[,2]),]
  result$medianGraph <- medianGraph
  return(result)
}
