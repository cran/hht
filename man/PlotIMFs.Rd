\name{PlotIMFs}
\alias{PlotIMFs}
\title{Display IMFs}
\description{
This function displays IMFs generated using \code{\link{Sig2IMF}}, \code{\link{EEMDCompile}.} or \code{\link{EEMDResift}}}
\usage{
PlotIMFs(sig, time.span = NULL, imf.list = NULL, original.signal = TRUE, 
    residue = TRUE, fit.line = FALSE, lwd = 1, cex = 1, ...)
}
\arguments{
\item{sig}{Data structure returned by \code{\link{Sig2IMF}}, \code{\link{EEMDCompile}}, or \code{\link{EEMDResift}}.}
\item{time.span}{Time span over which to plot IMFs.  \code{NULL} will draw the entire signal.}
\item{imf.list}{Which IMFs to plot, \code{NULL} plots all of them.}
\item{original.signal}{whether or not to plot the original signal.}
\item{residue}{whether to plot the residue of the EMD method.}
\item{fit.line}{whether to add a red line to the original signal trace showing how much of the original signal is contained in the selected IMFs and/or residual.}
\item{lwd}{Line weight.}
\item{cex}{Text size.}
\item{...}{Pass additional graphics parameters to IMF plotter}
}

\details{
This function plots the IMF decomposition of a signal.
It can show the original signal and also the residue left over when the IMFs are removed from the signal.
The plotter can use data from both EMD and EEMD runs.
When it plots EEMD data, it shows the averaged IMFs from the trials processed by \code{\link{EEMDCompile}}.
}

\note{
It is very important to inspect the IMF set prior to rendering Hilbert spectrograms.
Oftentimes, problems with the EMD are obvious when the IMFs are plotted.
The \code{fit.line} option can help with this.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{HHGramImage}}}
\examples{
data(PortFosterEvent)

#Run EMD
emd.result <- Sig2IMF(sig, tt, sm = "polynomial")

#Plot the first 4 IMFs of the EEMD of a signal.
time.span <- c(5, 10)
imf.list <- 1:4
original.signal <- TRUE
residue <- TRUE

PlotIMFs(emd.result, time.span, imf.list, original.signal, residue)

#Check how much contribution IMFs 2 and 3 make to the complete signal
imf.list <- c(2, 3)
fit.line <- TRUE
PlotIMFs(emd.result, time.span, imf.list, original.signal, residue, fit.line)
}

\keyword{ts}
\keyword{hplot}
