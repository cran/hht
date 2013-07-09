\name{EvolutiveFFT}
\alias{EvolutiveFFT}
\title{
Calculate the evolutive Fourier spectrogram.
}
\description{
Generates the evolutive Fourier spectrogram of a signal, and returns it for use in \code{\link{FTGramImage}}.
}
\usage{
EvolutiveFFT(sig, dt, ft, freq.span, taper = 0.05)
}
\arguments{
  \item{sig}{
  Signal to analyze.
}
  \item{dt}{
  Sample rate (must be constant).
}
  \item{ft}{
  Fourier transform input parameters
  \itemize{
      \item{ft$nfft}{  The frequency resolution, should be in powers of 2}
      \item{ft$ns}{  Number of samples in a window}
      \item{ft$nov}{  Number of samples to overlap, must be less than \code{ft$ns}}}
}
  \item{freq.span}{
  Frequency range to return.
}
  \item{taper}{
  Amount of cosine taper to apply.
}
}
\details{
This is an internal function and users will likely not call it directly.
}
\value{
   \item{z}{Power spectrum}
   \item{y}{Frequency}
   \item{x}{Time}
   \item{original.signal}{The input signal}
   \item{tt}{Sample times based on input sample rate \code{dt}}
}
\references{
Jonathan M. Lees (2012). RSEIS: Seismic Time Series Analysis Tools. R
package version 3.1-3.
}
\author{
Daniel C. Bowman \email{daniel.bowman@unc.edu}, Jonathan M. Lees
}
\note{
This is a modification of the \code{evolfft} function in the \code{RSEIS} package.
}

\keyword{ts}
