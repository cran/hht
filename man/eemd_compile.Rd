\name{eemd_compile}
\alias{eemd_compile}
\title{Process EEMD results}
\description{
This function compiles individual trial files from an EEMD run, allowing other functions to plot IMFs and Hilbert spectrograms of the data.
}
\usage{
eemd_compile(trials_dir, trials, nimf)  
}
\arguments{
\item{trials_dir}{Directory where previously generated EEMD trial files are stored.}
\item{trials}{Number of trial files to read.  This will warn users if the number of requested trials is greater than the number of files in the directory.}
\item{nimf}{Number of IMFs to record, IMFs past this number will not be saved.}
}
\details{
The EEMD algorithm can generate hundreds of files, resulting in massive amounts of data.
The \code{eemd_compile} function processes these files, generating an averaged IMF set and compiling the Hilbert spectrogram of each EMD run.
The output of \code{eemd_compile} can be used in \code{\link{plot_imfs}} and \code{\link{hhspec_image}}.
The averaged IMF set from \code{eemd_compile} can be resifted using \code{\link{eemd_resift}}. 
}
\value{ 
\item{eemd_result}{The averaged IMF set and individual Hilbert spectra of EMD trials generated through EEMD.}
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{eemd}}
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
