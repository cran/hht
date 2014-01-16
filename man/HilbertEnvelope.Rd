\name{HilbertEnvelope}
\alias{HilbertEnvelope}
\title{
Instantaneous amplitude
}
\description{
Generates the instantaneous amplitude of an analytic signal given by \code{\link{HilbertTransform}}
}
\usage{
HilbertEnvelope(asig)
}
\arguments{
  \item{asig}{
The analytic signal returned by \code{\link{HilbertTransform}}
}
}
\value{
\item{envelope}{Instantaneous amplitude}
}
\author{
Daniel C. Bowman \email{daniel.bowman@unc.edu}
}

\seealso{
\code{\link{HilbertTransform}}, \code{\link{InstantaneousFrequency}}
}
\examples{
tt <- seq(1000) * 0.01
sig <- sin(4 * pi * tt) + sin(3.4 * pi * tt)
asig <- HilbertTransform(sig)
env <- HilbertEnvelope(asig)
plot(tt, sig, type = "l")
lines(tt, env, col = "red")
lines(tt, -env, col = "red")

}
\keyword{ts}
