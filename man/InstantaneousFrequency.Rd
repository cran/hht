\name{InstantaneousFrequency}
\alias{InstantaneousFrequency}
\title{
Derive instantaneous frequency
}
\description{
Calculates instantaneous frequency from an analytic signal.
}
\usage{
InstantaneousFrequency(asig, tt, method = "arctan", lag = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{asig}{
Analytic signal produced by \code{\link{HilbertTransform}}
}
  \item{tt}{
Sample times
}
  \item{method}{
How the instantaneous frequency is calculated. 
\code{"arctan"} uses the arctangent of the real and imaginary parts of the Hilbert transform, taking the numerical derivative of phase for frequency.
\code{"chain"} uses the analytical derivative of the arctangent function prior to performing the numerical calculation.
}
  \item{lag}{
  Differentiation lag, see the \code{diff} function in the \code{base} package.
}
}
\value{
\item{instfreq}{Instataneous frequency in 1/time}
}
\author{
Daniel C. Bowman \email{daniel.bowman@unc.edu}
}
\note{

The \code{"arctan"} method was adapted from the \code{hilbertspec} function in the \code{EMD} package.

!!IMPORTANT!!
The numeric differentiation may be unstable for certain signals.
For example, high frequency sinusoids near the Nyquist frequency can give inaccurate results when using the \code{"chain"} method.
When in doubt, use the \code{\link{PrecisionTester}} function to check your results!
}

\seealso{
\code{\link{PrecisionTester}}
}
\keyword{ts}
