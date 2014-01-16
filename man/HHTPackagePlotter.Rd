\name{HHTPackagePlotter}
\alias{HHTPackagePlotter}
\title{
Set up spectrogram figure
}
\description{
Sets up the figure window for \code{\link{HHGramImage}} and \code{\link{FTGramImage}}.
This is an internal function and will likely never be called by a user
}

\usage{
HHTPackagePlotter(img, trace, amp.span, img.x.lab, img.y.lab, fit.line = NULL, 
window = NULL, colormap = NULL, backcol = c(0, 0, 0), pretty = FALSE, 
grid = TRUE, colorbar = TRUE, opts = list())
}
\arguments{
  \item{img}{
  Fourier or Hilbert spectrogram image.
}
  \item{trace}{
  Time series corresponding to the spectrogram.
}

  \item{amp.span}{
  Amplitudes over which to plot.
}
  \item{img.x.lab}{
  Specifies the X axis label on the image part of the figure, defaults to "time"
}
  \item{img.y.lab}{
  Specifies the Y axis label on the image part of the figure, defaults to "frequency"
}
  \item{fit.line}{
  Plots a line corresponding to the IMF sum on the trace, if requested
}
  \item{window}{
  The Fourier window length, if applicable
}
  \item{colormap}{
  The image color map
}
  \item{backcol}{
  The background color of the image (what shows up for pixels with value \code{NA})
}
  \item{pretty}{
  Adjusts image axes to have nice values, see the \code{pretty} function in the \code{base} package included in R
}
  \item{grid}{
  Determines whether to plot grid lines on the spectrogram
}
  \item{colorbar}{
  Whether to plot a color bar for amplitude values
}
  \item{opts}{
  Other possible options passed from \code{\link{HHGramImage}} and \code{\link{FTGramImage}}}
  
}
\value{
INTERNAL
}
\author{
Daniel Bowman \email{daniel.bowman@unc.edu}
}
\keyword{hplot}
