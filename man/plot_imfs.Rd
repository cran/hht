\name{plot_imfs}
\alias{plot_imfs}
\title{Display IMFs}
\description{
This function displays IMFs generated using \code{\link{sig2imf}}, \code{\link{eemd_compile}.} or \code{\link{eemd_resift}}}
\usage{
plot_imfs(sig, time_span, imf_list, original_signal, 
    residue, fit_line=FALSE, lwd=1, cex=1)
}
\arguments{
\item{sig}{Data structure returned by \code{\link{sig2imf}}, \code{\link{eemd_compile}.} or \code{\link{eemd_resift}}.}
\item{time_span}{Time span over which to plot IMFs.  \code{c(0, -1)} will draw the entire signal.}
\item{imf_list}{Which IMFs to plot.}
\item{original_signal}{(boolean) whether or not to plot the original signal.}
\item{residue}{(boolean) whether to plot the residue of the EMD method.}
\item{fit_line}{(boolean) whether to add a red line to the original signal trace showing how much of the original signal is contained in the selected IMFs and/or residual.}
\item{lwd}{Line weight.}
\item{cex}{Text size.}
}

\details{
This function plots the IMF decomposition of a signal.
It can show the original signal and also the residue left over when the IMFs are removed from the signal.
The plotter can use data from both EMD and EEMD runs.
When it plots EEMD data, it shows the averaged IMFs from the trials processed by \code{\link{eemd_compile}}.
}

\note{
It is very important to inspect the IMF set prior to rendering Hilbert spectrograms.
Oftentimes, problems with the EMD are obvious when the IMFs are plotted.
The \code{fit_line} option can help with this.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{hhspec_image}}}
\examples{
data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=5
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20

#Run EMD
emd_result=sig2imf(sig, dt, emd_config)

#Plot the first 4 IMFs of the EEMD of a signal.
time_span=c(5, 10)
imf_list=1:4
original_signal=TRUE
residue=TRUE

plot_imfs(emd_result, time_span, imf_list, original_signal, residue)

#Check how much contribution IMFs 2 and 3 make to the complete signal.
imf_list=c(2, 3)
fit_line=TRUE
plot_imfs(emd_result, time_span, imf_list, original_signal, residue, fit_line)
}

\keyword{nonparametric}
