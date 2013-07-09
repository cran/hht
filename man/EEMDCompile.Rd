\name{EEMDCompile}
\alias{EEMDCompile}
\title{Process EEMD results}
\description{
This function compiles individual trial files from an EEMD run, allowing other functions to plot IMFs and Hilbert spectrograms of the data.
}
\usage{
EEMDCompile(trials.dir, trials, nimf)  
}
\arguments{
\item{trials.dir}{Directory where previously generated EEMD trial files are stored.}
\item{trials}{Number of trial files to read.  This will warn users if the number of requested trials is greater than the number of files in the directory.}
\item{nimf}{Number of IMFs per EMD run. IMFs past this number will not be saved.}
}
\details{
The EEMD algorithm can generate hundreds of files, resulting in massive amounts of data.
The \code{EEMDCompile} function processes these files, generating an averaged IMF set and compiling the Hilbert spectrogram of each EMD run.
The output of \code{EEMDCompile} can be used in \code{\link{PlotIMFs}} and \code{\link{HHGramImage}}.
The averaged IMF set from \code{EEMDCompile} can be resifted using \code{\link{EEMDResift}}. 
}
\value{ 
\item{EEMD.result}{The averaged IMF set and individual Hilbert spectra of EMD trials generated through EEMD.}
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{EEMD}}, \code{\link{CombineTrials}}
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
