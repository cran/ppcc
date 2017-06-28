##  File ppcc.R
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
#' @title Plotting Point Positions
#' @description Calculates plotting point positions according
#' to different authors
#' @param n numeric, the sample size
#' @param method a character string naming a valid method (see Details)
#' @return a vector of class numeric that contains the plotting positions
#' @details
#' The following methods can by selected:
#' 
#' \code{"Gringorton"} the plotting point positions are
#' calculated as
#' \deqn{m_i = \left(i - 0.44\right) / \left(n + 0.12\right)}{%
#' m[i] = (i - 0.44) / (n + 0.12)}
#'
#' \code{"Cunane"} the plotting point positions are
#' calculated as
#' \deqn{m_i = \left(i - 0.4\right) / \left(n + 0.2\right)}{%
#' m[i] = (i - 0.4) / (n + 0.2)}
#'
#' \code{"Blom"} the plotting point positions are
#' calculated as
#' \deqn{m_i = \left(i - 0.3175\right) / \left(n + 0.25\right)}{%
#' m[i] = (i - 0.3175) / (n + 0.25)}
#'
#' \code{"Filliben"} the order statistic medians are calculated as:
#' \deqn{
#' m_i = \left\{
#'       \begin{array}{l l}
#'       1 - 0.5^{1/n} & i = 1 \\
#'       \left(i - 0.3175\right)/\left(n + 0.365\right) & i = 2,\ldots, n - 1 \\
#'       0.5^{1/n} & i = n \\
#'       \end{array}\right.
#' }{%
#' m[1] = 1 - 0.5^(1/n); m[i] = (i - 0.3175) / (n + 0.365),
#' for i = 2,\ldots, n - 1; m[n] = 0.5^(1/n)
#' }
#'
#' \code{"ppoints"} R core's default plotting point positions
#' are calculated (see \code{\link{ppoints}}).
#' 
#' @importFrom stats ppoints
#' @export ppPositions
#' @keywords misc

ppPositions <- function(n, method= c("Gringorton", "Cunane",
                                   "Filliben", "Blom", "Weibull", "ppoints"))
{
    method <- match.arg(method)

    if (method == "Gringorton"){
        m <- numeric(n)
        for (i in 1:n){
            m[i] <- (i - 0.44) / (n + 0.12)
        }
    } else if (method == "Cunane") {
        m <- numeric(n)
        for (i in 1:n){
            m[i] <- (i - 0.4) / (n + 0.2)
        }
    } else if (method == "Filliben"){
        m <- numeric(n)
        m[1] <- 1 - 0.5^(1/n)
        for (i in 2:(n-1)){
            m[i] <- (i - 0.3175) / (n + 0.365)
        }
        m[n] <- 0.5^(1/n)

    } else if (method == "Blom") {
        m <- numeric(n)
        for (i in 1:n){
            m[i] <- (i - 0.375) / (n + 0.25)
        }

    } else if (method == "Weibull"){
        m <- numeric(n)
        for (i in 1:n){
            m[i] <- i / (n + 1)
        }

    #} # else if (method == "Nguyen") {
      #  n <- length(x)
      #  m <- numeric(n)
      #  gamma <- skewness(x)
      #  stopifnot(gamma >= -3 & gamma <= 3 & n >= 5 & n <= 100)
      #  for (i in 1:n){
      #      m[i] <- (i - 0.42) / (n + 0.3 * gamma + 0.05)
      #  }
   # }
    } else if (method == "ppoints"){
        m <- ppoints(n)
    }
       
    return(m)
}
