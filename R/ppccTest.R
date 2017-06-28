##  File ppccTest.R
##  Part of the R package ppcc
##
##  Copyright (C) 2017 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/
#
#' @title Probability Plot Correlation Coefficient Test
#' 
#' @description Performs the Probability Plot Correlation Coeffient
#' Test of Goodness-of-Fit
#' 
#' @concept Goodness-of-Fit
#' @concept Composite Hypothesis Test
#' @concept Homogeinity test
#' @param x a numeric vector of data values; NA values will be silently ignored.
#' @param qfn a character vector naming a valid quantile function
#' @param shape numeric, the shape parameter for the relevant distribution,
#' if applicable; defaults to NULL 
#' @param mc numeric, the number of Monte-Carlo replications, defaults to 10000
#' @param ppos character, the method for estimating plotting point positions,
#' default's to \code{NULL}, see Details for corresponding
#' defaults and \code{\link{ppPositions}} for available methods
#' @param \ldots further arguments, currently ignored
#' @return a list with class 'htest'
#' 
#' @details
#' Filliben (1975) suggested a probability plot correlation
#' coeffient test to test a sample for normality. The ppcc is defined as
#' the product moment correlation coefficient between the
#' ordered data \eqn{x_{(i)}}{x[(i)]} and the order statistic medians \eqn{M_{i}}{M[i]},
#' \deqn{
#' r =  \frac{\sum_{i = 1}^n \left(x_{(i)} - \bar{x} \right)~ \left(M_i - \bar{M}\right)}
#' {\sqrt{\sum_{i=1}^n \left(x_{(i)} - \bar{x}\right)^2 ~ \sum_{j = 1}^n \left(M_j - \bar{M} \right)^2}},}{%
#' r = Cov(X, M) / [Var(X) * Var(M)]^(1/2),}
#'
#' whereas the ordered statistic medians are related to the quantile function
#' of the standard normal distribution, \eqn{M_{i} = \phi^{-1} (m_i)}{M[i] = phi^(-1) (m[i])}.
#' The values of \eqn{m_i}{m[i]} are estimated by plotting-point position procedures
#' (see \code{\link{ppPositions}}).
#'
#' In this function the test is performed by Monte-Carlo simulation:
#' \enumerate{
#' \item Calculate quantile-quantile \eqn{\hat{r}}{r*} for the ordered
#' sample data \code{x}
#' and the specified \code{qfn} distribution (with \code{shape}, if applicable)
#' and given \code{ppos}.
#' \item Draw \code{n} (pseudo) random deviates from the specified \code{qfn}
#' distribution, where \code{n} is the sample size of \code{x}.
#' \item Calculate quantile-quantile \eqn{r_i}{r[i]} for the random
#' deviates and the specified \code{qfn} distribution with given \code{ppos}.
#' \item Repeat step 2 and 3 for \eqn{i = \left\{1, 2, \ldots, mc\right\}}{i = 1, 2, ..., mc}.
#' \item Calculate \eqn{S = \sum_{i=1}^n \mathrm{sgn}(\hat{r} - r_i)}{S = \sum(sgn(r* - r[i]))} with
#' sgn the sign-function.
#' \item The estimated \eqn{p}-value is \eqn{p = S / mc}{p = S / mc}.
#' }
#' 
#' The probability plot correlation coeffient is invariant for \code{location}
#' and \code{scale}. Therefore, the null hypothesis is a
#' composite hypothesis, e.g. \eqn{H0: X \in N(\mu, \sigma),
#' ~~ \mu \in R,~~ \sigma \in R_{>0}}{%
#' H0: X ~ N(\mu, \sigma), \mu = [-\infty, \infty], \sigma > 0}.
#' Furthermore, distributions with one (additional) specified
#' \code{shape} parameter can be tested.
#' 
#' The magnitude of \eqn{\hat{r}}{r*} depends on the selected method for
#' plotting-point positions (see \code{\link{ppPositions}})
#' and the sample size. Several authors extended Filliben's method to
#' assess the goodness-of-fit to other distributions, whereas theoretical
#' quantiles were used as opposed to Filliben's medians.
#'
#' The default plotting positions (see \code{\link{ppPositions}})
#' depend on the selected \code{qfn}.
#' 
#' Distributions with none or one single scale parameter
#' that can be tested:
#' \tabular{llll}{
#' Argument \tab Function \tab Default pppos \tab Reference \cr
#' \code{qunif} \tab Uniform \tab \code{Weibull} \tab Vogel and Kroll (1989) \cr
#' \code{qexp} \tab Exponential \tab \code{Gringorton} \tab \cr
#' \code{qgumbel} \tab Gumbel \tab \code{Gringorton} \tab Vogel (1986) \cr
#' \code{qrayleigh} \tab Rayleigh \tab \code{Gringorton} \tab \cr
#' }
#'
#' Distributions with \code{location} and \code{scale} parameters
#' that can be tested:
#' \tabular{llll}{
#' Argument \tab Function \tab Default pppos \tab Reference \cr
#' \code{qnorm} \tab Normal \tab \code{Blom} \tab Looney and Gulledge (1985) \cr
#' \code{qlnorm} \tab log-Normal \tab \code{Blom} \tab Vogel and Kroll (1989)\cr
#' \code{qcauchy} \tab Cauchy \tab \code{Gringorton} \tab \cr
#' \code{qlogis} \tab Logistic \tab \code{Blom} \tab \cr
#' }
#'
#' If Blom's plotting position is used for \code{qnorm}, than the ppcc-test
#' is related to the Shapiro-Francia
#' normality test (Royston 1993), where \eqn{W' = r^2}. See
#' \code{\link[nortest]{sf.test}} and \code{example(ppccTest)}.
#'
#' Distributions with additional \code{shape} parameters
#' that can be tested:
#' \tabular{llll}{
#' Argument \tab Function \tab Default pppos \tab Reference \cr
#' \code{qweibull} \tab Weibull \tab \code{Gringorton} \tab \cr
#' \code{qpearson3} \tab Pearson III \tab \code{Blom} \tab Vogel and McMartin (1991)\cr
#' \code{qgev} \tab GEV \tab \code{Cunane} \tab Chowdhury et al. (1991) \cr
#' \code{qkappa2} \tab two-param. Kappa Dist. \tab \code{Gringorton} \tab  \cr
#' \code{qglogis} \tab Generalized Logistic \tab \code{Gringorton} \tab  \cr
#' }
#'
#' If \code{qfn = qpearson3} and \code{shape = 0} is selected, the
#' \code{qnorm} distribution is used. If \code{qfn = qgev} and
#' \code{shape = 0}, the \code{qgumbel} distribution is used.
#' If \code{qfn = qglogis} and \code{shape = 0} is selected, the
#' \code{qglogis} distribution is used.
#' 
#' @note
#' As the \code{pvalue} is estimated through a Monte-Carlo simulation,
#' the results depend on the selected seed (see \code{\link{set.seed}})
#' and the total number of replicates (\code{mc}).
#'
#' The default of \code{mc = 10000} re-runs is sufficient for
#' testing the composite hypothesis on levels of \eqn{\alpha = [0.1, 0.05]}.
#' If a level of \eqn{\alpha = 0.01} is desired, than larger sizes
#' of re-runs (e.g. \code{mc = 100000}) might be required. 
#'
#' @seealso
#' \code{\link{qqplot}}, \code{\link{qqnorm}}, \code{\link{ppoints}},
#' \code{\link{ppPositions}}, \code{\link{Normal}},
#' \code{\link{Lognormal}}, \code{\link{Uniform}}, \code{\link{Exponential}}, \code{\link{Cauchy}},
#' \code{\link{Logistic}}, \code{\link[VGAM]{qgumbel}}, \code{\link{Weibull}},
#' \code{\link[VGAM]{qgev}}.
#' 
#' @references
#' J. U. Chowdhury, J. R. Stedinger, L.-H. Lu (1991), Goodness-of-Fit Tests
#' for Regional Generalized Extreme Value Flood Distributions,
#' \emph{Water Resources Research} 27, 1765--1776.
#' 
#' J. J. Filliben (1975), The Probability Plot Correlation
#' Coefficient Test for Normality, \emph{Technometrics} 17, 111--117.
#'
#' S. Kim, H. Shin, T. Kim, J.-H. Heo (2010), Derivation of the Probability
#' Plot Correlation Coefficient Test Statistics for the Generalized Logistic
#' Distribution. Intern. Workshop Adv. in Stat. Hydrol., May 23 - 25, 2010
#' Taormina.
#' 
#' S. W. Looney, T. R. Gulledge (1985), Use of Correlation Coefficient
#' with Normal Probability Plots, \emph{The American Statistician} 39,
#' 75--79.
#'
#' P. W. Mielke (1973), Another family of distributions for
#' describing and analyzing precipitation data.
#' \emph{Journal of Applied Meteorology} 12, 275--280.
#'
#' P. Royston, P. (1993), A pocket-calculator algorithm for the
#' Shapiro-Francia test for non-normality: an application to
#' medicine. \emph{Statistics in Medicine} 12, 181-184.
#' 
#' R. M. Vogel (1986), The Probability Plot Correlation Coefficient
#' Test for the Normal, Lognormal, and Gumbel Distributional
#' Hypotheses, \emph{Water Resources Research} 22, 587--590.
#' 
#' R. M. Vogel, C. N. Kroll (1989), Low-flow frequency analysis using
#' probability-plot correlation coefficients, \emph{Journal of Water
#' Resources Planning and Management} 115, 338--357.
#'
#' R. M. Vogel, D. E. McMartin (1991), Probability Plot Goodness-of-Fit
#' and Skewness Estimation Procedures for the Pearson Type 3
#' Distribution, \emph{Water Resources Research} 27, 3149--3158.
#' 
#' @examples
#' ## Filliben (1975, p.116)
#' ## Note: Filliben's result was 0.98538
#' ## decimal accuracy in 1975 is assumed to be less than in 2017
#' x <- c(6, 1, -4, 8, -2, 5, 0)
#' set.seed(100)
#' ppccTest(x, "qnorm", ppos="Filliben")
#' ## p between .75 and .9
#' ## see Table 1 of Filliben (1975, p.113)
#' ##
#' set.seed(100)
#' ## Note: default plotting position for
#' ## qnorm is ppos ="Blom"
#' ppccTest(x, "qnorm")
#' ## p between .75 and .9
#' ## see Table 2 of Looney and Gulledge (1985, p.78)
#' ##
#' ## 
#' set.seed(300)
#' x <- rnorm(30)
#' qn <- ppccTest(x, "qnorm")
#' qn
#' ## p between .5 and .75
#' ## see Table 2 for n = 30 of Looney and Gulledge (1985, p.78)
#' ##
#' ## Compare with Shapiro-Francia test
#' if(require(nortest)){
#'    sn <- sf.test(x)
#'    print(sn)
#'    W <- sn$statistic
#'    rr <- qn$statistic^2
#'    names(W) <- NULL
#'    names(rr) <- NULL
#'    print(all.equal(W, rr))
#' }
#' ppccTest(x, "qunif")
#' ppccTest(x, "qlnorm")
#' old <- par()
#' par(mfrow=c(1,3))
#' xlab <- "Theoretical Quantiles"
#' ylab <- "Empirical Quantiles"
#' qqplot(x = qnorm(ppPositions(30, "Blom")),
#'        y = x, xlab=xlab, ylab=ylab, main = "Normal q-q-plot")
#' qqplot(x = qunif(ppPositions(30, "Weibull")),
#'        y = x, xlab=xlab, ylab=ylab, main = "Uniform q-q-plot")
#' qqplot(x = qlnorm(ppPositions(30, "Blom")),
#'        y = x, xlab=xlab, ylab=ylab, main = "log-Normal q-q-plot")
#' par(old)
#' ##
#' if (require(VGAM)){
#' set.seed(300)
#' x <- rgumbel(30)
#' gu <- ppccTest(x, "qgumbel")
#' print(gu)
#' 1000 * (1 -  gu$statistic)
#' }
#' ##
#' ## see Table 2 for n = 30 of Vogel (1986, p.589) 
#' ## for n = 30 and Si = 0.5, the critical value is 16.9 
#' ##
#' set.seed(200)
#' x <- runif(30)
#' un <- ppccTest(x, "qunif")
#' print(un)
#' 1000 * (1 - un$statistic)
#' ##
#' ## see Table 1 for n = 30 of Vogel and Kroll (1989, p.343)
#' ## for n = 30 and Si = 0.5, the critical value is 10.5
#' ##
#' set.seed(200)
#' x <- rweibull(30, shape = 2.5)
#' ppccTest(x, "qweibull", shape=2.5)
#' ppccTest(x, "qweibull", shape=1.5)
#' ##
#' if (require(VGAM)){
#' set.seed(200)
#' x <- rgev(30, shape = -0.2)
#' ev <- ppccTest(x, "qgev", shape=-0.2)
#' print(ev)
#' 1000 * (1 - ev$statistic)
#' ##
#' ## see Table 3 for n = 30 and shape = -0.2
#' ## of Chowdhury et al. (1991, p.1770)
#' ## The tabulated critical value is 80.
#' }
#' 
#' @keywords htest
#' @keywords npar
#' 
#' @useDynLib 'ppcc', .registration = TRUE
#'
#' @export ppccTest
#'
#' @importFrom stats qnorm qweibull qgamma qexp qunif qlnorm qlogis
#' @importFrom stats qcauchy
###################################################################
ppccTest <-
    function(x, qfn =  c("qnorm", "qlnorm", "qunif", "qexp",
                          "qcauchy", "qlogis", "qgumbel",
                          "qweibull", "qpearson3", "qgev",
                          "qkappa2", "qrayleigh", "qglogis"),
             shape = NULL, ppos = NULL, mc = 10000, ...)
{
   
    qfn <- match.arg(qfn)
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    if(n < 1L)
        stop("not enough 'x' data")
    if(!is.character(qfn))
        stop("'qfn' must be a string naming a valid function")

    if (qfn == "qnorm"){
        ##
        if (is.null(ppos)) ppos <- "Blom" #Default, user may change to Filliben
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from a Normal distribution")

        ## Get ppcc first
        ## theoretical quantiles
        q <- qnorm(pe)
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n, res = res)$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
        PVALUE <- .C("ppcctest_norm", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval)$pval

     ##   if (ppos == "Blom"){
            ## From Eq. 5a of Heo et al. (2008, p.4)
     ##       f <- 1 / exp(1.29 + 0.283 * 0.05 +
     ##           (0.887 - 0.751 * 0.05 + 3.21 * 0.05 * 0.05) * log(n)) 
     ##       ESTIMATE <- 1 - f
     ##       names(ESTIMATE) <- "r(alpha < 0.05)"
     ##   } else {
            ESTIMATE <- NULL
     ##   }
    } else  if (qfn == "qlnorm"){
        if (is.null(ppos)) ppos <- "Blom"  ## Blom's plotting points
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from a log-Normal distribution")

        ## Get ppcc first
        ## theoretical quantiles
        q <- qlnorm(pe)
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)

        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
        
        PVALUE <- .C("ppcctest_lnorm", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval )$pval
        ESTIMATE <- NULL
    } else if (qfn == "qunif"){
        ## Weibull plotting positions
        ## R. M. Vogel, C. N. Kroll (1989), Low-flow frequency analysis
        ## using probability-plot correlation coefficients,
        ## Journal of Water Resources Planning and Management 115, 338--357.
        ## see p. 342
        if (is.null(ppos)) ppos <- "Weibull"  ## Weibulls plotting points
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from a Uniform distribution")

        ## Get ppcc first
        ## theoretical quantiles
        q <- qunif(pe)
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
 
        PVALUE <- .C("ppcctest_unif", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval )$pval
        ESTIMATE <- NULL
    
    } else if (qfn == "qrayleigh"){
        ## Gringorton's plotting positions
        if (is.null(ppos)) ppos <- "Gringorton"  ## 
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from a Rayleigh distribution")

        ## Get ppcc first
        ## theoretical quantiles
        q <- sqrt(-2 * log(1 - pe))
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
 
        PVALUE <- .C("ppcctest_rayleigh", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval )$pval
        ESTIMATE <- NULL
    
    }   else if (qfn == "qexp"){
        ## Exponential distribution
        if (is.null(ppos)) ppos <- "Gringorton"  ## Weibulls plotting points
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from an Exponential distribution")

        ## Get ppcc first
        ## theoretical quantiles
        q <- qexp(pe)
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
 
        PVALUE <- .C("ppcctest_exp", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval )$pval
        ESTIMATE <- NULL

    } else if (qfn == "qlogis"){
        ## Logistic distribution
        if (is.null(ppos)) ppos <- "Blom" # similar to Normal
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from a Logistic distribution")

        ## Get ppcc first
        ## theoretical quantiles
        q <- qlogis(pe)
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
 
        PVALUE <- .C("ppcctest_logis", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval )$pval
        ESTIMATE <- NULL
        
    } else if (qfn == "qcauchy"){
        ## Cauchy distribution
        if (is.null(ppos)) ppos <- "Gringorton" 
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from a Cauchy distribution")

        ## Get ppcc first
        ## theoretical quantiles
        q <- qcauchy(pe)
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
 
        PVALUE <- .C("ppcctest_cauchy", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval )$pval
              ESTIMATE <- NULL

    } else if (qfn == "qgumbel"){
        ## Gumbel distribution
        if (is.null(ppos)) ppos <- "Gringorton"
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME, " differs from a Gumbel distribution")

        ## Get ppcc first
        ## theoretical quantiles for Gumbel
        q <- -log(-log(pe))
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
 
        PVALUE <- .C("ppcctest_gumbel", p = pe, ppcc = r, n = n, nmc = mc,
                     pval = pval )$pval
        ESTIMATE <- NULL
        
   }  else if (qfn == "qpearson3"){
        ## pearson3 distribution

        if(is.null(shape)){
            stop("If 'qpearson3' is selected, the parameter
                  'shape' must be specified.")
        }
        
        if (is.null(ppos)) ppos <- "Blom"
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME,
                       " differs from a Pearson-3 distribution\nwith shape ",
                        shape)

        ## Get ppcc first
        ## theoretical quantiles for the Pearson-3  distribution
        #q <- y(pe, ...)
        if (shape == 0){
            q <- qnorm(pe, 0, 1)
        } else if (shape > 0){
            alpha <- 4 / shape^2
            beta <- 1 / 2 * abs(shape)
            q <- -alpha * beta * qgamma(pe, alpha, scale = beta)
        } else {
            alpha <- 4 / shape^2
            beta <- 1 / 2 * abs(shape)
            q <- alpha * beta - qgamma(1 - pe, alpha, scale = beta)
        }
        q[pe == 0 & shape > 0] <- -2 / shape
        q[pe == 1 & shape < 0] <- -2 / shape
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
        shape <- as.double(shape)
        if (shape == 0){
            PVALUE <- .C("ppcctest_norm", p = pe, ppcc = r, n = n, nmc = mc,
                         pval = pval )$pval
        } else {
            PVALUE <- .C("ppcctest_pearson3", p = pe, ppcc = r, shape = shape,
                         sn = n, nmc = mc, pval = pval )$pval
        }
        ESTIMATE <- NULL
        
    }  else if (qfn == "qweibull"){
        ## Weibull distribution

        if(is.null(shape)){
            stop("If 'qweibull' is selected, the parameter
                  'shape' must be specified.")
        }
        
        if (is.null(ppos)) ppos <- "Gringorton"
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME,
                       " differs from a Weibull distribution\nwith shape ",
                        shape)

        ## Get ppcc first
        ## theoretical quantiles for the Weibull distribution
        #q <- y(pe, ...)
        q <- qweibull(pe, shape=shape, scale=1)
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
        shape <- as.double(shape)
        PVALUE <- .C("ppcctest_weibull", p = pe, ppcc = r, shape = shape,
                     sn = n, nmc = mc, pval = pval )$pval
        ESTIMATE <- NULL
        
    } else if (qfn == "qgev"){
        ## GEV distribution

        if(is.null(shape)){
            stop("If 'qgev' is selected, the parameter
                  'shape' must be specified.")
        }
        
        if (is.null(ppos)) ppos <- "Cunane"
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME,
                       " differs from a Generalized Extreme Value distribution\nwith shape ",
                        shape)

        ## Get ppcc first
        ## theoretical quantiles for the GEV distribution
        #q <- y(pe, ...)
        if (shape == 0){
             ## theoretical quantiles for Gumbel
            q <- -log(-log(pe))
        } else {
            q <- 1 / shape * (1 - ( -log(pe))^shape)
        }
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
        shape <- as.double(shape)
        if (shape == 0){
            PVALUE <- .C("ppcctest_gumbel", p = pe, ppcc = r,
                         sn = n, nmc = mc, pval = pval )$pval
        } else {
            PVALUE <- .C("ppcctest_gev", p = pe, ppcc = r, shape = shape,
                         sn = n, nmc = mc, pval = pval )$pval
        }
        ESTIMATE <- NULL


    } else if (qfn == "qglogis"){
        ## GEV distribution

        if(is.null(shape)){
            stop("If 'qglogis' is selected, the parameter
                  'shape' must be specified.")
        }
        
        if (is.null(ppos)) ppos <- "Gringorton"
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME,
                       " differs from a Generalized Logistic distribution\nwith shape ",
                        shape)

        ## Get ppcc first
        ## theoretical quantiles for the Generalized Logistic distribution
        #q <- y(pe, ...)
        if (shape == 0){
             ## theoretical quantiles
            q <- qlogis(pe)
        } else {
            q <- (1 - exp(-shape * log(pe / (1 - pe)))) / shape
        }
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
        shape <- as.double(shape)
        if (shape == 0){
            PVALUE <- .C("ppcctest_logis", p = pe, ppcc = r,
                         sn = n, nmc = mc, pval = pval )$pval
        } else {
            PVALUE <- .C("ppcctest_glogis", p = pe, ppcc = r, shape = shape,
                         sn = n, nmc = mc, pval = pval )$pval
        }
        ESTIMATE <- NULL

    } else if (qfn == "qkappa2"){
        ## Mielke's Kappa distribution

        if(is.null(shape)){
            stop("If 'qkappa2' is selected, the parameter
                  'shape' must be specified.")
        }
        
        if (is.null(ppos)) ppos <- "Gringorton"
        pe <- ppPositions(n, ppos)
        ALTERNATIVE <- paste0(DNAME,
                       " differs from Mielke's Kappa distribution\nwith shape ",
                        shape)

        ## Get ppcc first
        ## theoretical quantiles for the Kappa distribution
        #q <- y(pe, ...)
        if (shape == 0){
            stop("Kappa distribution for 'shape = 0' is not defined")
        } else {
            q <- (shape * pe^shape / (1 - pe^shape))^( 1 / shape)
        }
        ## empirical quantiles
        qe <- sort(x)
        
        qe <- as.double(qe)
        q <- as.double(q)
        n <- as.integer(n)
        res <- as.double(1.0)
        
        STATISTIC <- .C("pmcor", x = qe, y = q, n = n,
                        res = res )$res
        r <- as.double(STATISTIC)
        mc <- as.integer(mc)
        pval <- as.double(1.0)
        pe <- as.double(pe)
        shape <- as.double(shape)
        
        PVALUE <- .C("ppcctest_kappa2", p = pe, ppcc = r, shape = shape,
                     sn = n, nmc = mc, pval = pval )$pval     
        ESTIMATE <- NULL
    }
    
    
    names(STATISTIC) <- "ppcc"
    PARAMETER <- n
    names(PARAMETER) <- "n"
    METHOD <- "Probability Plot Correlation Coefficient Test"
    ans <- list(data.name = DNAME, method = METHOD, p.value = PVALUE,
                statistic = STATISTIC, alternative = ALTERNATIVE,
                estimate = ESTIMATE, parameter = PARAMETER)
    class(ans) <- "htest"  
    return(ans)
}
