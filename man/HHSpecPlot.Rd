\name{HHSpecPlot}
\alias{HHSpecPlot}
\title{Display Hilbert periodogram}
\description{
This function displays the Hilbert periodogram, with options to plot individual IMFs and also the Fourier periodogram for comparison.
}
\usage{
HHSpecPlot(hspec, freq.span = NULL, scaling = "none", imf.list = NULL, 
    show.total = TRUE, show.fourier = FALSE, scale.fourier = FALSE, 
    show.imfs = FALSE, legend = TRUE, ...)
}
\arguments{
\item{hspec}{Data structure returned by \code{\link{HHSpectrum}}}
\item{freq.span}{Frequency range to plot, \code{NULL} plots all of them}
\item{scaling}{Amplitude scaling, can be \code{"log"} (log 10), \code{"sqrt"} (square root), defaults to \code{"none"}.}
\item{imf.list}{Which IMFs to plot, requires \code{show.imfs = TRUE}.}
\item{show.total}{Show the ensemble Hilbert spectrogram}
\item{show.fourier}{Show the Fourier periodogram}
\item{scale.fourier}{Scale Fourier and Hilbert spectra to each other for easier comparison}
\item{show.imfs}{Plot individual IMF spectra}
\item{legend}{Determines whether or not a legend is shown}
\item{\dots}{This function supports some optional parameters as well:
\itemize{
    \item{xlab}{ - X axis label}
    \item{ylab}{ - Y axis label}
    \item{legend.location}{ - where to put the legend}
    \item{total.col}{ - color of ensemble Hilbert periodogram}
    \item{total.lwd}{ - lwd of ensemble Hilbert periodogram}
    \item{total.lty}{ - lty of ensemble Hilbert periodogram}
    \item{imf.cols}{ - colors of IMF periodogram}
    \item{imf.lwd}{ - lwds of IMF periodogram}
    \item{imf.lty}{ - ltys of IMF periodogram}
    \item{fourier.col}{ - color of Fourier periodogram}
    \item{fourier.lwd}{ - lwd of Fourier periodogram}
    \item{fourier.lty}{ - lty of Fourier periodogram}
    \item{main}{ - figure title}}
}}

\details{
This function plots the Hilbert periodogram of a signal, with options to show periodograms of individual IMFs.
You can also plot a simple Fourier periodogram for comparison.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{HHSpectrum}}, \code{\link{HHGramImage}}}
\examples{

#Here we see how the EMD produces a dyadic filter bank for uniform random noise
#The frequency distributions of all but the first IMF display a Chi-Square distribution
#See Huang, N. E. & Wu, Z. 
#A review on Hilbert-Huang Transform: Method and its applications to geophysical studies.
#Reviews of Geophysics, 2008, 46, RG2006

#The EMD of this signal may take a couple of minutes to run

set.seed(628)
sig  <-  runif(10000)
tt  <-  seq_len(length(sig)) * 0.01

\dontrun{emd.result  <-  Sig2IMF(sig, tt)}

dfreq  <-  0.1
\dontrun{hspec  <-  HHSpectrum(emd.result, dfreq)}

\dontrun{HHSpecPlot(hspec, show.imfs = TRUE, 
imf.list = 1:10, show.total = TRUE, scaling = "sqrt", 
imf.lwd = rep(2, 10), total.lty = 3)}
}

\keyword{ts}
\keyword{hplot}
