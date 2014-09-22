\name{PrecisionTester}
\alias{PrecisionTester}
\title{
Test numerically determined instantaneous frequency against exact instantaneous frequency
}
\description{
This function compares the performance of \code{\link{InstantaneousFrequency}} against signals of known instantaneous frequency.
The known signal is of the form
\deqn{ x(t) = a\sin(\omega_{1} + \varphi_{1}) + b\sin(\omega_{2} + \varphi_{2}) + c}{a * \sin(omega_1 + phi_1) + b * sin(omega_2 + phi_2) + c}
One can create quite complicated signals by choosing the various amplitude, frequency, and phase constants.
}
\usage{
PrecisionTester(tt = seq(0, 10, by = 0.01), method = "arctan", lag = 1, 
    a = 1, b = 1, c = 1, omega.1 = 2 * pi, omega.2 = 4 * pi, 
    phi.1 = 0, phi.2 = pi/6, plot.signal = TRUE, 
    plot.instfreq = TRUE, plot.error = TRUE, new.device = TRUE, ...)
}
\arguments{
  \item{tt}{
Sample times.
}
  \item{method}{
How the numeric instantaneous frequency is calculated, see \code{\link{InstantaneousFrequency}}
}
  \item{lag}{
Differentiation lag, see the \code{diff} function in the \code{base} package.
}
  \item{a}{
Amplitude coefficient for the first sinusoid.
}
  \item{b}{
Amplitude coefficient for the second sinusoid.
}
  \item{c}{
DC shift
}
  \item{omega.1}{
Frequency of the first sinusoid.
}
  \item{omega.2}{
Frequency of the second sinusoid.
}
  \item{phi.1}{
Phase shift of the first sinusoid.
}
  \item{phi.2}{
Phase shift of the second sinusoid.
}
  \item{plot.signal}{
Whether to show the time series.
}
  \item{plot.instfreq}{
Whether to show the instantaneous frequencies, comparing the numerical and analytical result.
}
  \item{plot.error}{
Whether to show the difference between the numerical and analytical result.
}
  \item{new.device}{
  Whether to open each plot as a new plot window (defaults to \code{TRUE}).  However, Sweave doesn't like \code{dev.new()}.
  If you want to use \code{PrecisionTester} in Sweave, be sure that \code{new.device = FALSE}}
  \item{\dots}{
Plotting parameters.
}
}
\details{
For a description of how the exact analytical frequency is derived, see \url{http://www.unc.edu/\%7Ehaksaeng/hht/analytic_instantaneous_freq.pdf}
}
\value{
\item{instfreq$sig}{The time series}
\item{instfreq$analytic}{The exact instantaneous frequency}
\item{instfreq$numeric}{The numerically-derived instantaneous frequency from \code{\link{InstantaneousFrequency}}}
}
\author{
Daniel C. Bowman \email{daniel.bowman@unc.edu}
}

\examples{

#Simple signal
tt <- seq(0, 10, by = 0.01)
a <- 1
b <- 0
c <- 0
omega.1 <- 30 * pi
omega.2 <- 0
phi.1 <- 0
phi.2 <- 0
PrecisionTester(tt, method = "arctan", lag = 1, a, b, c, 
    omega.1, omega.2, phi.1, phi.2)

#That was nice - what happens if we use the "chain" method...?

PrecisionTester(tt, method = "chain", lag = 1, a, b, c, 
    omega.1, omega.2, phi.1, phi.2)

#Big problems!  Let's increase the sample rate

tt <- seq(0, 10, by = 0.0005)
PrecisionTester(tt, method = "chain", lag = 1, a, b, c, 
omega.1, omega.2, phi.1, phi.2)

#That's better

#Frequency modulations caused by signal that is not symmetric about 0

tt <- seq(0, 10, by = 0.01)
a <- 1
b <- 0
c <- 0.25
omega.1 <- 2 * pi
omega.2 <- 0
phi.1 <- 0
phi.2 <- 0

PrecisionTester(tt, method = "arctan", lag = 1, a, b, c, 
omega.1, omega.2, phi.1, phi.2)

#Non-uniform sample rate
set.seed(628)
tt <- sort(runif(500, 0, 10))
a <- 1
b <- 0
c <- 0
omega.1 <- 2 * pi
omega.2 <- 0
phi.1 <- 0
phi.2 <- 0

PrecisionTester(tt, method = "arctan", lag = 1, a, b, c, 
omega.1, omega.2, phi.1, phi.2)
}

\seealso{\code{\link{InstantaneousFrequency}}}

\keyword{ts}
