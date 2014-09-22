\name{HilbertTransform}
\alias{HilbertTransform}
\title{
The Hilbert transform
}
\description{
Creates the analytic signal using the Hilbert transform.
}
\usage{
HilbertTransform(sig)
}
\arguments{
  \item{sig}{
Signal to transform.
}
}
\details{
Creates the real and imaginary parts of a signal.
}
\value{
\item{asig}{Analytic signal}
}
\author{
Daniel C. Bowman \email{daniel.bowman@unc.edu}
}

\seealso{
\code{\link{HilbertEnvelope}}, \code{\link{InstantaneousFrequency}}
}
\examples{
tt <- seq(1000) * 0.01
sig <- sin(pi * tt)
asig <- HilbertTransform(sig)
plot(tt, sig, xlim = c(0, 12))
lines(tt, Re(asig), col = "green")
lines(tt, Im(asig), col = "red")
legend("topright", col = c("black", "green", "red"), 
lty = c(NA, 1, 1), pch = c(1, NA, NA), 
legend = c("Signal", "Real", "Imaginary"))

}
\keyword{ts}
