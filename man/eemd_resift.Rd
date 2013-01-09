\name{eemd_resift}
\alias{eemd_resift}
\title{Resift averaged IMFs from EEMD}
\description{
Averaged IMFs produced by EEMD may not satisfy the strict definition of an IMF, and therefore they may not have meaningful Hilbert spectrograms.
Huang and Wu (2008) suggest another round of sifting to ensure that the averaged IMFs are made to satisfy the IMF definition.
This function resifts the averaged IMF set and saves the results based on rules described in the input \code{resift_rule}.
}
\usage{
eemd_resift(eemd_result, emd_config, resift_rule)  
}
\arguments{
\item{eemd_result}{The averaged IMF set and individual Hilbert spectra of EMD trials generated through EEMD.}
\item{emd_config}{Configuration parameters for the resifting emd algorithm, see the \code{emd} function in the \code{EEMD} package for a detailed description of each option:}
\itemize{
    \item{\code{emd_config$max_sift}}{ maximum number of IMF sifts}
    \item{\code{emd_config$max_imf}}{ maximum number of IMFs that can be returned}
    \item{\code{emd_config$tol}}{ sifting stop criterion}
    \item{\code{emd_config$stop_rule}}{ EMD stop rules}
    \item{\code{emd_config$boundary}}{ how the start and stop of the time series are handled in the splining process}
    \item{\code{emd_config$sm}}{ spline smoothing}
    \item{\code{emd_config$spar}}{ smoothing parameter}
    \item{\code{emd_config$weight}}{ spline weight}
}
\item{resift_rule}{How the resifting algorithm chooses which IMF to save}
\itemize{
    \item{Integer - Which IMF in the resifted set will be saved (so if resift_rule=1, the first IMF will be saved, the rest will be discarded)}
    \item{``last'' - The last IMF will be saved (not terribly useful)}
    \item{``max_var'' - The IMF with the most variance will be saved.  This will get the most ``significant'' IMF out of each resifted set.}
    \item{``all'' - Every single new IMF generated from resifting the averaged IMFs will be saved.  There may be a lot of them!}}
}
\details{
The function \code{\link{eemd_compile}} generates a list of averaged IMFs from EEMD trials.
These averaged IMFs often do not satisfy the definition of an IMF, usually because some of them are mixtures of different time scales.
This is a consequence of the noise perturbation method of EEMD, but it complicates the attempt to create a meaningful Hilbert spectrogram from the averaged IMF set.
The resifting algorithm takes each averaged IMF and performs EMD, thereby splitting each one into multiple ``sub-IMFs'', each of which satisfy the strict definition of an IMF.
The question then is: which of these sub-IMFs best represent the averaged IMF
The most rigorous solution is to set \code{resift_rule} to \code{"all"}, but that tends to make a large number of sub-IMFs, many with very low amplitude.
Another solution is to accept the sub-IMF with the most variance, as that probably represents the fundamental information content of the original averaged IMF.
}
\value{ 
\item{resift_result}{The resifted results of the averaged IMF set and the individual Hilbert spectra of each resifted IMF.}
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{eemd}}, \code{\link{eemd_compile}}
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

trials=10
nimf=10
noise_amp=6.4e-07
trials_dir="test"

set.seed(628)

#Run EEMD (this may take some time)
\dontrun{eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)}

#Compile the results
\dontrun{eemd_result=eemd_compile(trials_dir, trials, nimf)}


resift_rule="max_var"
\dontrun{resift_result=eemd_resift(eemd_result, emd_config, resift_rule)}

#Plot the IMFs
time_span=c(5, 10)
imf_list=1:3
os=TRUE
res=TRUE
\dontrun{plot_imfs(resift_result, time_span, imf_list, os, res)}
}
\keyword{nonparametric}
