#' Simulate & test mediation methods as in Barfield et al. (2017)
#'
#' Simulate & test mediation methods as in Barfield et al. (2017). If multiple genes are simulated, only 1st
#' is tested. For each simulation, count proportion of p-values < \code{alpha}, and at end calculate mean proportion.
#'
#' @param med.fnm Quoted name of mediation function to test. The function must accept parameters and output matrix
#' with mediation p-value in column \code{"EMY.p"}, like \code{\link{hitman}}.
#' @param b1t2.v Numeric vector of values that both theta2 (\code{"t2"}) and beta1 (\code{"b1"}) take.
#' @param nsamp Number of samples.
#' @param ngene Number of genes other than that of primary interest to simulate.
#' @param ... Sent to \code{med.fcn}
#' @inheritParams ezlimma:::sim_fisher
#' @return Matrix with proportion of significant calls for every combination of \code{t2} and \code{b1}.
#' @references Barfield R, Shen J, Just AC, Vokonas PS, Schwartz J, Baccarelli AA, VanderWeele TJ, Lin X.
#' Testing for the indirect effect under the null for genome-wide mediation analyses. Genet Epidemiol.
#' 2017 Dec;41(8):824-833.
#' @export

sim_barfield <- function(med.fnm, b1t2.v=c(0, 0.14, 0.39), alpha=0.05, nsamp=50, nsim=10**4, ngene=0,
                              seed=0, verbose=TRUE, ...){
  med.fcn <- eval(parse(text=med.fnm))
  stopifnot(ngene == 0 || med.fnm == "hitman")

  #t = theta; b = beta
  t0 <- t1 <- t3 <- b0 <- b2 <- 0.14
  prop.sig.arr <- array(NA, dim=c(length(b1t2.v), length(b1t2.v), nsim),
                        dimnames=list(paste0("t2_", b1t2.v), paste0("b1_", b1t2.v), paste0("sim_", 1:nsim)))
  # sim <- 1
  # t2 <- b1 <- 0.4
  # calculate E(M) & E(Y)
  # simulate
  set.seed(seed)
  for (sim in 1:nsim){
    for (t2 in b1t2.v){
      for (b1 in b1t2.v){
        x <- stats::rnorm(n=nsamp)
        a <- stats::rnorm(n=nsamp)

        # eq 1; E(M1)
        em1 <- b0 + b1*a + b2*x
        # normal here ~ log-normal or inv chi sq on counts
        m1 <- stats::rnorm(n=nsamp, mean=em1)

        # eq 3; E(Y)
        ey <- t0 + t1*a + t2*m1 + t3*x
        y <- stats::rnorm(n=nsamp, mean=ey)
        names(y) <- paste0("s", 1:length(y))

        if (ngene >= 1){
          em.other <- b0+b2*x
          m.other.mat <- matrix(stats::rnorm(n=nsamp*ngene, mean=em.other, sd=2), nrow=ngene, ncol=nsamp, byrow = TRUE)
          med.mat <- rbind(m1, m.other.mat)
          dimnames(med.mat) <- list(paste0("m", 1:nrow(med.mat)), paste0("s", 1:ncol(med.mat)))
          med.res <- hitman(E=a, M=med.mat, Y=y, covariates = x, verbose=FALSE)
          prop.sig.arr[paste0("t2_", t2), paste0("b1_", b1), paste0("sim_", sim)] <- med.res["m1", "EMY.p"] < alpha
        } else {
          # ngene = 0 || not hitman
          if (med.fnm == "lotman"){
            med.res <- lotman(E=a, M=m1, Y=y, covariates = x, verbose = FALSE)
          } else {
            med.res <- med.fcn(E=a, M=m1, Y=y, covariates = x, ...)
          }
          prop.sig.arr[paste0("t2_", t2), paste0("b1_", b1), paste0("sim_", sim)] <- med.res[1, "EMY.p"] < alpha
        }
      }
    }
    if (verbose && sim %% 100 == 0) message("sim: ", sim)
  }
  # summarize
  prop.sig.mat <- apply(prop.sig.arr, MARGIN=c(1,2), FUN=mean)
}
