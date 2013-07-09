\name{CombineTrials}
\alias{CombineTrials}
\title{Gather EEMD trial files}
\description{
This function gathers trial files from multiple directories, renumbers them, and saves them to a single directory for processing using \code{\link{EEMDCompile}}.
}
\usage{
CombineTrials(in.dirs, out.dir, copy=TRUE)  
}
\arguments{
\item{in.dirs}{Directories containing trial file sets from one EEMD run.}
\item{out.dir}{Directory in which to save all trial files.}
\item{copy}{Copy files (\code{TRUE}) or move them (\code{FALSE}).}
}
\details{
Parallel processing is an efficient method for running EEMD.
However, this will result in several directories, each with trial files numbered from 1 to N.
These files cannot simply be copied together into the same directory, because then they would overwrite each other.
This function gathers all trial files in multiple directories, renumbers them, and saves them in a different directory.
}
\value{ 
The trial files are saved in the directory specified by \code{out.dir}.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{EEMD}}, \code{\link{EEMDCompile}}
}
\examples{
#Suppose you have run 3 different EEMD sets of 100 trials each 
#and saved the results in EEMD1, EEMD2, EEMD3, respectively:
in.dirs <- c("/home/user/EEMD1", "/home/user/EEMD2/", "/home/user/EEMD3")
out.dir <- "/home/user/all.trials"
\dontrun{CombineTrials(in.dirs, out.dir)}
#Now all your trials should be located in /home/user/all.trials, 
#numbered 1 through 300
}
\keyword{ts}
