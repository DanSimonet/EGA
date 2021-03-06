#'  From the EGA's structure to CFA.
#'
#' \code{CFA} Verifies the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @param ega.obj An EGA object.
#' @param estimator The estimator used in the confirmatory factor analysis. 'WLSMV' is the estimator of choice for ordinal variables. 'ML' or 'WLS' for interval variables.
#'
#' @param plot.CFA Logical. Should the CFA structure with its standardized loadings be plot?
#' @param data A dataframe with the variables to be used in the analysis.
#' @param ... Arguments passed to ’cfa’ in lavaan.
#' @author Hudson F. Golino <hfgolino at gmail.com>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24])
#' cfa.wmt <- CFA(ega.obj = ega.wmt, estimator = 'WLSMV', plot.CFA = TRUE, data = wmt2)
#'
#' ega.intel <- EGA(data = intelligenceBattery[,8:66])
#' cfa.intel <- CFA(ega.obj = ega.intel, estimator = 'WLSMV', plot.CFA = TRUE, data = intelligenceBattery[,8:66])
#'
#' \dontrun{
#' CFA(a, estimator = 'WLSMV', data = data, ...)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap.
#' @export


CFA <- function(ega.obj, estimator, plot.CFA = TRUE, data, layout = "spring", ...) {
  if(!require(lavaan)) {
    message("installing the 'lavaan' package")
    install.packages("lavaan")
  }

  if(!require(semPlot)) {
    message("installing the 'semPlot' package")
    install.packages("semPlot")
  }

    strct <- split(ega.obj$dim.variables[, 1], list(ega.obj$dim.variables[, 2]))
    names(strct) <- paste("Fat", labels(strct))
    model.ega <- paste(names(strct), " =~ ", lapply(strct, function(x) paste(print(x), collapse = " + ")), collapse = " \n ")
    fit.mod.ega <- cfa(model = model.ega, estimator = estimator, orthogonal = FALSE, se = "standard", test = "satorra-bentler",
        data = data, ...)
    summary.cfa <- summary(fit.mod.ega, fit.measures = TRUE)
    fit.measures.cfa <- fitMeasures(fit.mod.ega, fit.measures = c("chisq", "df", "pvalue", "cfi", "rmsea", "gfi", "nfi"))

    if (plot.CFA == TRUE) {
        plot.cfa <- semPaths(fit.mod.ega, title = FALSE, label.cex = 0.5, sizeLat = 5, sizeMan = 3, edge.label.cex = 0.4, minimum = 0.1,
            sizeInt = 0.5, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
            "std", cut = 0.5)
    }

    cfa <- list()
    cfa$fit <- fit.mod.ega
    cfa$summary <- summary.cfa
    cfa$fit.measures <- fit.measures.cfa
    return(cfa)
}
