\name{eemd}
\alias{eemd}
\title{Ensemble Empirical Mode Decomposition}
\description{
This function performs ensemble empirical mode decomposition (EEMD).
}
\usage{
eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir=NULL)  
}
\arguments{
\item{sig}{a time series to be decomposed (vector)}
\item{dt}{The sample rate of \code{sig}}
\item{trials}{Number of times to run EMD}
\item{nimf}{Number of IMFs to record, IMFs past this number will not be saved}
\item{noise_amp}{Amplitude of white noise to use in denoising algorithm}
\item{emd_config}{Configuration parameters for EMD algorithm, see \code{emd} function in the \code{EMD} package for a detailed description of what each option does.}
\itemize{
\item{\code{emd_config$tol}}{ sifting stop criterion}
\item{\code{emd_config$stop_rule}}{ EMD stop rules}
\item{\code{emd_config$boundary}}{ how the start and stop of the time series are handled in the splining process}
\item{\code{emd_config$sm}}{ spline smoothing}
\item{\code{emd_config$spar}}{ smoothing parameter}
\item{\code{emd_config$weight}}{ spline weight}}
\item{trials_dir}{Directory where EEMD trial files will be stored, defaults to ``trials.''  This will create a directory if none exists.}
}
\details{
This function performs ensemble empirical mode decomposition, a noise assisted version of the EMD algorithm.
The EEMD works by adding a certain amplitude of white noise to a time series, decomposing it via EMD, and saving the result.
If this is done enough times, the averages of the noise perturbed IMFs will approach the ``true'' IMF set.
The EEMD can ameliorate mode mixing and intermittency problems (see references section).

This EEMD algorithm creates a directory \code{trials_dir} and saves each EMD trial into this directory.
The number of trials is defined using \code{trials}.
The trial files in this directory can then be processed using \code{\link{eemd_compile}} to produce the averaged IMF set, or to plot the Hilbert spectrogram of the data.
Keep in mind that the EEMD is an expensive algorithm and may take significant time to run.
}
\value{ 
\item{emd_result}{The result of each individual EMD trial. This is saved directly to file in directory \code{trials_dir} (i.e. it is not returned by \code{eemd}.)}
}
\references{
Wu, Z. A. and Huang, N. E. (2009) Ensemble empirical mode decomposition: {A} noise assisted data analysis method.
\emph{Advances in Adaptive Data Analysis}, \bold{1}, 1-41.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{sig2imf}}, \code{\link{eemd_compile}}, \code{\link{plot_imfs}}.
}
\examples{
data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

trials=10
nimf=10
noise_amp=6.4e-07
trials_dir="test"

set.seed(628)
#Run EEMD (this may take some time)
\dontrun{eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)}

#Compile the results
\dontrun{eemd_result=eemd_compile(trials_dir, trials, nimf)}

#Plot the IMFs
time_span=c(5, 10)
imf_list=1:3
os=TRUE
res=TRUE
\dontrun{plot_imfs(eemd_result, time_span, imf_list, os, res)}
}
\keyword{nonparametric}
