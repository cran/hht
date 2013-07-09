\name{EEMDResift}
\alias{EEMDResift}
\title{Resift averaged IMFs from EEMD}
\description{
Averaged IMFs produced by EEMD may not satisfy the strict definition of an IMF, and therefore they may not have meaningful Hilbert spectrograms.
Huang and Wu (2008) suggest another round of sifting to ensure that the averaged IMFs are made to satisfy the IMF definition.
This function resifts the averaged IMF set and saves the results based on rules described in the input \code{resift.rule}.
}
\usage{
EEMDResift(EEMD.result, resift.rule, spectral.method = "arctan", 
    diff.lag = 1, tol = 5, max.sift = 200, stop.rule = "type5", 
    boundary = "wave", sm = "none", smlevels = c(1), 
    spar = NULL, max.imf = 100, interm = NULL)  
}
\arguments{
\item{EEMD.result}{The averaged IMF set and individual Hilbert spectra of EMD trials generated through EEMD.}
\item{resift.rule}{How the resifting algorithm chooses which IMF to save
\itemize{
    \item{Integer - Which IMF in the resifted set will be saved (so if resift.rule=1, the first IMF will be saved, the rest will be discarded)}
    \item{``last'' - The last IMF will be saved (not terribly useful)}
    \item{``max.var'' - The IMF with the most variance will be saved.  This will get the most ``significant'' IMF out of each resifted set.}
    \item{``all'' - Every single new IMF generated from resifting the averaged IMFs will be saved.  There may be a lot of them!}}}
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
}

\details{
The function \code{\link{EEMDCompile}} generates a list of averaged IMFs from EEMD trials.
These averaged IMFs often do not satisfy the definition of an IMF, usually because some of them are mixtures of different time scales.
This is a consequence of the noise perturbation method of EEMD, but it complicates the attempt to create a meaningful Hilbert spectrogram from the averaged IMF set.
The resifting algorithm takes each averaged IMF and performs EMD, thereby splitting each one into multiple ``sub-IMFs'', each of which satisfy the strict definition of an IMF.
The question then is: which of these sub-IMFs best represent the averaged IMF?
The most rigorous solution is to set \code{resift.rule} to \code{"all"}, but that tends to make a large number of sub-IMFs, many with very low amplitude.
Another solution is to accept the sub-IMF with the most variance, as that probably represents the fundamental information content of the original averaged IMF.
}
\value{ 
\item{resift.result}{The resifted results of the averaged IMF set and the individual Hilbert spectra of each resifted IMF.}
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{EEMD}}, \code{\link{EEMDCompile}}
}
\examples{
data(PortFosterEvent)

trials=10
nimf=10
noise.amp=6.4e-07
trials.dir="test"

set.seed(628)

#Run EEMD (this may take some time)
\dontrun{EEMD(sig, tt, noise.amp, trials, nimf, noise.amp, trials.dir = trials.dir)}

#Compile the results
\dontrun{EEMD.result <- EEMDCompile(trials.dir, trials, nimf)}


resift.rule="max.var"
\dontrun{resift.result <- EEMDResift(EEMD.result, resift.rule)}

#Plot the IMFs
time.span=c(5, 10)
imf.list=1:3
os=TRUE
res=TRUE
\dontrun{PlotIMFs(resift.result, time.span, imf.list, os, res)}
}
\keyword{ts}
