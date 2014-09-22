\name{CEEMD}
\alias{CEEMD}
\title{Complete Ensemble Empirical Mode Decomposition}
\description{
This function implements the complete ensemble empirical mode decomposition (CEEMD) algorithm.
}
\usage{
CEEMD(sig, tt, noise.amp, trials, verbose = TRUE, 
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
    Each row in the array corresponds to the noise series added for that particular trial during the CEEMD run.
    By default, \code{noise.array = NULL}.}
}
\details{
This function performs the complete ensemble empirical mode decomposition, a noise assisted empirical mode decomposition algorithm.
The CEEMD works by adding a certain amplitude of white noise to a time series, decomposing it via EMD, and saving the result.
In contrast to the Ensemble Empirical Mode Decomposition (EEMD) method, the CEEMD also ensures that the IMF set is quasi-complete and orthogonal.
The CEEMD can ameliorate mode mixing and intermittency problems.
Keep in mind that the CEEMD is a computationally expensive algorithm and may take significant time to run.
}
\value{ 
\item{ceemd.result}{The final result of the CEEMD algorithm}.
}
\references{
Torres, M. E., Colominas, M. A., Schlotthauer, G., Flandrin, P. (2011).  A complete ensemble empirical mode decomposition with adaptive noise.
\emph{2011 IEEE International Conference on Acoustics, Speech, and Signal Processing}, pp.4144-4147, doi: 10.1109/ICASSP.2011.5947265.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{EEMD}}, \code{\link{Sig2IMF}}, \code{\link{PlotIMFs}}.
}
\examples{

\dontrun{

data(PortFosterEvent)
noise.amp <- 6.4e-07
trials <- 100

ceemd.result <- CEEMD(sig, tt, noise.amp, trials)
PlotIMFs(ceemd.result, imf.list = 1:6, time.span = c(5, 10))
}

}
\keyword{ts}
