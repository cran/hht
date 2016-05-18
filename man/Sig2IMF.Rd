\name{Sig2IMF}
\alias{Sig2IMF}
\title{Empirical Mode Decomposition wrapper}
\description{
This function wraps the \code{emd} function in the \code{EMD} package.
\code{Sig2IMF} is used in \code{\link{EEMD}} and others.
}
\usage{
Sig2IMF(sig, tt, spectral.method = "arctan", diff.lag = 1, stop.rule = "type5", 
    tol = 5, boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, 
    max.sift = 200, max.imf = 100, interm = NULL)
}

\arguments{
\item{sig}{a time series to be decomposed (vector)}
\item{tt}{A vector of sample times for \code{sig}}
\item{spectral.method}{defines how to calculate instantaneous frequency - whether to use the arctangent of the analytic signal with numeric differentiation (``arctan'')
    or the result of the chain rule applied to the arctangent, then numerically differentiated ("chain");  see \code{\link{InstantaneousFrequency}}.}
\item{diff.lag}{specifies if you want to do naive differentiation (\code{diff.lag} = 1), central difference method (\code{diff.lag = 2} or higher difference methods (\code{diff.lag} > 2) 
     to determine instantaneous frequency; see \code{\link{InstantaneousFrequency}}.}
\item{stop.rule}{As quoted from the EMD package documentation:  ''The stop rule of sifting. 
     The type1 stop rule indicates that absolute values 
     of envelope mean must be less than the user-specified tolerance level in the sense
     that the local average of upper and lower envelope is zero. The stopping rules
     type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
     of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
     and Wu (2008), respectively.''}
\item{tol}{Determines what value is used to stop the sifting - this will depend on which stop rule you use.}
\item{boundary}{how the beginning and end of the signal are handled}
\item{sm}{Specifies how the signal envelope is constructed, see Kim et al, 2012.}
\item{smlevels}{Specifies what level of the IMF is obtained by smoothing other than interpolation, see EMD package documentation}
\item{spar}{User-defined smoothing parameter for spline, kernel, or local polynomial smoothing.}
\item{max.sift}{How many sifts are allowed - if this value is exceeded the IMF is returned as-is.}
\item{max.imf}{Maximum number of IMFs allowed.}
\item{interm}{Specifies vector of periods to be excluded from IMFs to cope with mode mixing.}
}
\details{
This function configures and performs empirical mode decomposition using the \code{emd} function in the \code{EMD} package.
}
\value{ 
\item{emd.result}{The intrinsic mode functions (IMFs), instantaneous frequencies, and instantaneous amplitudes of \code{sig}.}
}
\references{
Kim, D., Kim, K. and Oh, H.-S. (2012) Extending the scope of empirical mode decomposition by smoothing.
\emph{EURASIP Journal on Advances in Signal Processing}, \bold{2012}, 168.

Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.

Huang, N. E. and Wu Z. A. (2008) A review on Hilbert-Huang Transform: Method and its applications to geophysical studies.
\emph{Reviews of Geophysics}, \bold{46}, RG2006.
}

\seealso{
\code{\link{EEMD}}, \code{\link{PlotIMFs}}
}
\examples{
data(PortFosterEvent)

#Run EMD
emd.result=Sig2IMF(sig, tt)

#Display IMFs

time.span <- c(5, 10)
imf.list <- 1:3
original.signal <- TRUE
residue <- TRUE

PlotIMFs(emd.result, time.span, imf.list, original.signal, residue)

#Plot spectrogram
sdt <- tt[2] - tt[1]
dfreq <- 0.25
freq.span <- c(0, 25)
hgram <- HHRender(emd.result, sdt, dfreq, freq.span = freq.span, verbose = FALSE)

time.span <- c(4, 10)
freq.span <- c(0, 25)
amp.span <- c(0.000001, 0.00001)
HHGramImage(hgram, time.span = time.span, 
freq.span = freq.span, amp.span = amp.span,
pretty = TRUE)
}
\keyword{ts}
