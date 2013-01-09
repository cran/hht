\name{sig2imf}
\alias{sig2imf}
\title{Empirical Mode Decomposition wrapper}
\description{
This function wraps the \code{emd} function in the \code{EMD} package.
\code{sig2imf} is used in \code{\link{eemd}} and others.
}
\usage{
sig2imf(sig, dt, emd_config)
}
\arguments{
\item{sig}{a time series to be decomposed (vector)}
\item{dt}{The sample rate of \code{sig}}
\item{emd_config}{Configuration information for the EMD algorithm, see the \code{emd} function in the \code{EMD} package for a detailed description of what each option does.}
\itemize{
    \item{\code{emd_config$max_sift}}{ maximum number of IMF sifts}
    \item{\code{emd_config$max_imf}}{ maximum number of IMFs that can be returned}
    \item{\code{emd_config$tol}}{ sifting stop criterion}
    \item{\code{emd_config$stop_rule}}{ EMD stop rules}
    \item{\code{emd_config$boundary}}{ how the start and stop of the time series are handled in the splining process}
    \item{\code{emd_config$sm}}{ spline smoothing}
    \item{\code{emd_config$spar}}{ smoothing parameter}
    \item{\code{emd_config$weight}}{ spline weight}}
}
\details{
This function configures and performs empirical mode decomposition.
}
\value{ 
\item{emd_result}{A data structure for input into \code{\link{eemd_compile}}, \code{\link{plot_imfs}}, and \code{\link{hhtransform}}.}
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.

Huang, N. E. and Wu Z. A. (2008) A review on Hilbert-Huang Transform: Method and its applications to geophysical studies.
\emph{Reviews of Geophysics}, \bold{46}, RG2006.
}

\seealso{
\code{\link{eemd}}, \code{\link{plot_imfs}}
}
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

#Display IMFs

time_span=c(5, 10)
imf_list=1:3
original_signal=TRUE
residue=TRUE

plot_imfs(emd_result, time_span, imf_list, original_signal, residue)

#Get Hilbert transform
emd_result=hhtransform(emd_result)

#Plot spectrogram
max_freq=25
freq_step=0.05
hspec=hh_render(emd_result, max_freq, freq_step)

freq_span=c(0, 25)
amp_span=c(0.000001, 0.00001)
hhspec_image(hspec, time_span, freq_span, amp_span)
}
\keyword{nonparametric}
