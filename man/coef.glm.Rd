% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coefglm.R
\name{coef.glm}
\alias{coef.glm}
\title{Coefficient and Bias From Generalized Linear Model Fit}
\usage{
\method{coef}{glm}(object, type = c("raw", "bias", "corrected"), ...)
}
\arguments{
\item{object}{An object with similar structure as returned from 
\code{\link[stats:glm]{glm}}.}

\item{type}{A character scalar, being one of \code{'raw'}, 
\code{'bias'}, and \code{'corrected'}. \code{'raw'} returns the 
usual coefficients; \code{'bias'} returns the estimated leading order bias;
and \code{'corrected'} returns the coefficients after subtracting
the estimated bias.}

\item{...}{Not used.}
}
\value{
A numeric vector of requested components. 
 When \code{type='corrected'}, the \code{'bias'} attribute
 will be set to a numeric vector of estimated biases 
 being subtracted from the raw coefficients.
}
\description{
\code{coef.glm} is a \code{S3} method computing the raw coefficients, 
leading order bias, and bias corrected coefficients from a 
generalized linear model (glm) fit.
}
\examples{
 x=1:30
 y=rpois(30L, x/10)
 glmfit=glm(y~x, poisson('log'))
 coef(glmfit)
 coef(glmfit, type='bias')
 coef(glmfit, type='corrected')
 
}
\references{
McCullagh and Nelder (1989) "Generalized Linear Models", 2nd edition. London: Chapman and Hall.

Cordeiro and McCullah (1991) "Bias Correction in Generalized Linear Models", \emph{Journal of the Royal Statistical Society: Series B}, \bold{53}, 629-643.
}
\author{
Long Qu <rtistician@gmail.com>
}
\concept{
bias correction
}
\keyword{methods}
