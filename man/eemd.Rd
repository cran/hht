\name{EEMD}
\alias{EEMD}
\title{Ensemble Empirical Mode Decomposition}
\description{
This function performs ensemble empirical mode decomposition (EEMD).
}
\usage{
EEMD(sig, tt, noise.amp, trials, nimf, trials.dir = NULL, verbose = TRUE, 
    spectral.method = "arctan", diff.lag = 1, tol = 5, max.sift = 200,
    stop.rule = "type5", boundary = "wave", sm = "none",
    smlevels = c(1), spar = NULL, max.imf = 100, interm = NULL, 
    noise.type = "gaussian", noise.array = NULL)
}
\arguments{
\item{sig}{a time series to be decomposed (vector)}
\item{tt}{The sample times of \code{sig}}
\item{noise.amp}{Amplitude of white noise to use in denoising algorithm}
\item{trials}{Number of times to run EMD}
\item{nimf}{Number of IMFs to record, IMFs past this number will not be saved}
\item{trials.dir}{Directory where EEMD trial files will be stored, defaults to ``trials.''  This will create a directory if none exists.}
\item{verbose}{If TRUE, notify when each trial is complete}
\item{spectral.method}{See \code{\link{Sig2IMF}}.}
\item{diff.lag}{See \code{\link{Sig2IMF}}.}
\item{tol}{See \code{\link{Sig2IMF}}.}
\item{max.sift}{See \code{\link{Sig2IMF}}.}
\item{stop.rule}{See \code{\link{Sig2IMF}}.}
\item{boundary}{See \code{\link{Sig2IMF}}.}
\item{sm}{See \code{\link{Sig2IMF}}.}
\item{smlevels}{See \code{\link{Sig2IMF}}.}
\item{spar}{See \code{\link{Sig2IMF}}.}
\item{max.imf}{See \code{\link{Sig2IMF}}.}
\item{interm}{See \code{\link{Sig2IMF}}.}
\item{noise.type}{If unspecified or \code{gaussian}, produce a Gaussian noise series with length \code{length(sig)} and standard deviation \code{noise.amp}.
    If \code{uniform}, produce a uniform random distribution with length \code{length(sig)} and maximum absolute value of \code{noise.amp}.
    If \code{custom}, then use a custom noise array as defined in input parameter \code{noise.array} (see below).}
\item{noise.array}{If \code{noise.type = "custom"}, this array must be a TRIALS x LENGTH(TT) collection of time series to be used in the place of uniform or gaussian noise.
    Each row in the array corresponds to the noise series added for that particular trial during the EEMD run.
    By default, \code{noise.array = NULL}.}
}
\details{
This function performs ensemble empirical mode decomposition, a noise assisted version of the EMD algorithm.
The EEMD works by adding a certain amplitude of white noise to a time series, decomposing it via EMD, and saving the result.
If this is done enough times, the averages of the noise perturbed IMFs will approach the ``true'' IMF set.
The EEMD can ameliorate mode mixing and intermittency problems (see references section).

This EEMD algorithm creates a directory \code{trials.dir} and saves each EMD trial into this directory.
The number of trials is defined using \code{trials}.
The trial files in this directory can then be processed using \code{\link{EEMDCompile}} to produce the averaged IMF set, or to plot the Hilbert spectrogram of the data.
Keep in mind that the EEMD is an expensive algorithm and may take significant time to run.
}
\value{ 
\item{emd.result}{The result of each individual EMD trial. This is saved directly to files in directory \code{trials.dir} (i.e. it is not returned by \code{EEMD}.)}
}
\references{
Wu, Z. A. and Huang, N. E. (2009) Ensemble empirical mode decomposition: {A} noise assisted data analysis method.
\emph{Advances in Adaptive Data Analysis}, \bold{1}, 1-41.
}

\note{Previous versions of this function used a uniform random noise distribution (i.e. \code{runif}) to generate the noise time series.
The default noise time series is now Gaussian in accordance with existing EEMD literature.}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{Sig2IMF}}, \code{\link{CombineTrials}}, \code{\link{EEMDCompile}}, \code{\link{PlotIMFs}}.
}
\examples{
data(PortFosterEvent)
trials <- 10
nimf <- 10
noise.amp <- 6.4e-07
trials.dir <- "test"

set.seed(628)
#Run EEMD (this may take some time)
\dontrun{EEMD(sig, tt, noise.amp, trials, nimf, trials.dir = trials.dir)}

#Compile the results
\dontrun{EEMD.result <- EEMDCompile(trials.dir, trials, nimf)}

#Plot the IMFs
time.span <- c(5, 10)
imf.list <- 1:3
os <- TRUE
res <- TRUE
\dontrun{PlotIMFs(EEMD.result, time.span, imf.list, os, res)}
}
\keyword{ts}
