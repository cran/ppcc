% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppPositions.R
\name{ppPositions}
\alias{ppPositions}
\title{Plotting Point Positions}
\usage{
ppPositions(
  n,
  method = c("Gringorton", "Cunane", "Filliben", "Blom", "Weibull", "ppoints")
)
}
\arguments{
\item{n}{numeric, the sample size}

\item{method}{a character string naming a valid method (see Details)}
}
\value{
a vector of class numeric that contains the plotting positions
}
\description{
Calculates plotting point positions according
to different authors
}
\details{
The following methods can by selected:

\code{"Gringorton"} the plotting point positions are
calculated as
\deqn{m_i = \left(i - 0.44\right) / \left(n + 0.12\right)}{%
m[i] = (i - 0.44) / (n + 0.12)}

\code{"Cunane"} the plotting point positions are
calculated as
\deqn{m_i = \left(i - 0.4\right) / \left(n + 0.2\right)}{%
m[i] = (i - 0.4) / (n + 0.2)}

\code{"Blom"} the plotting point positions are
calculated as
\deqn{m_i = \left(i - 0.3175\right) / \left(n + 0.25\right)}{%
m[i] = (i - 0.3175) / (n + 0.25)}

\code{"Filliben"} the order statistic medians are calculated as:
\deqn{
m_i = \left\{
      \begin{array}{l l}
      1 - 0.5^{1/n} & i = 1 \\
      \left(i - 0.3175\right)/\left(n + 0.365\right) & i = 2,\ldots, n - 1 \\
      0.5^{1/n} & i = n \\
      \end{array}\right.
}{%
m[1] = 1 - 0.5^(1/n); m[i] = (i - 0.3175) / (n + 0.365),
for i = 2,\ldots, n - 1; m[n] = 0.5^(1/n)
}

\code{"ppoints"} R core's default plotting point positions
are calculated (see \code{\link{ppoints}}).
}
\keyword{misc}
