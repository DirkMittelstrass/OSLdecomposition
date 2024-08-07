% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_PhotoCrosssections.R
\name{plot_PhotoCrosssections}
\alias{plot_PhotoCrosssections}
\title{Plot comparison of CW-OSL component photoionisation cross sections of different models}
\usage{
plot_PhotoCrosssections(
  fit.list,
  stimulation.intensity = NULL,
  stimulation.wavelength = NULL,
  K.selected = NULL,
  title = NULL,
  hide.plot = FALSE,
  filename = NULL
)
}
\arguments{
\item{fit.list}{\link{list} (\strong{required}):
Output object of \link{fit_OSLcurve}. The object must be created with the setting \code{output.complex = TRUE}.}

\item{stimulation.intensity}{\link{numeric} (\emph{optional}):
Intensity of optical stimulation in \emph{mW / cm²}. Used to calculate the photoionisation cross sections.
If not given, the input value for \link{fit_OSLcurve} is used}

\item{stimulation.wavelength}{\link{numeric} (\emph{optional}):
Wavelength of optical stimulation in nm. Used to calculate the photoionisation cross sections.
If not given, the input value for \link{fit_OSLcurve} is used}

\item{K.selected}{\link{numeric} (\emph{optional}):
Draws a red rectangle around the \code{K = K.selected} row, thus highlighting the model of choice.}

\item{title}{\link{character} (\emph{with default}):
Plot title. Set \code{title = NULL} for no title.}

\item{hide.plot}{\link{logical} (\emph{with default}):
If true, plot is not drawn but can still be saved as file or caught by \code{A <- plot_PhotoCrosssections(...)}.
If caught, the plot can be drawn manually for example by using \link[gridExtra:arrangeGrob]{gridExtra::grid.arrange}.}

\item{filename}{\link{character} (\emph{optional}):
File name or path to save the plot as image. If just a file name is given, the image is
saved in the working directory. The image type is chosen by the file ending. Both, vector images
as well as pixel images are possible. Allowed are \code{.pdf}, \code{.eps}, \code{.svg} (vector graphics),
\code{.jpg}, \code{.png}, \code{.bmp} (pixel graphics) and more, see \link[ggplot2:ggsave]{ggplot2::ggsave}.}
}
\value{
An invisible \link[ggplot2:ggplot]{ggplot2::ggplot} object containing the diagram will returned. "Invisible" means, the no value
will be returned (e.g. no console printout) if the function is not assigned to a variable via \verb{<-}.
If the function is assigned, the returned object can be further manipulated with
\link[ggplot2:ggplot2-package]{ggplot2::ggplot2-package} methods or manually drawn by various functions like for example \link[gridExtra:arrangeGrob]{gridExtra::grid.arrange}.
}
\description{
This function takes the \code{output.complex = TRUE} output of \link{fit_OSLcurve} and draws the
photoionisation cross sections of different models in relation to each other.
If a stimulation wavelength between 465 and 480 nm was chosen,
the photoionisation cross sections are also set in relation to literature values
from Singarayer and Bailey (2003), Jain et al. (2003) and Durcan and Duller (2011).

The photoionisation cross section ranges of the reference components are defined as following:

\tabular{lll}{
\strong{Component} \tab \strong{Lower limit (cm^2)} \tab \strong{Upper limit (cm^2)}\cr
Ultrafast \tab 1e-16 \tab 1e-15 \cr
Fast \tab 1.9e-17 \tab 3.1e-17 \cr
Medium \tab 3e-18 \tab 9e-18 \cr
Slow1 \tab 1e-18 \tab 1.85e-18 \cr
Slow2 \tab 1.1e-19 \tab 4e-19 \cr
Slow3 \tab 1e-20 \tab 4.67e-20 \cr
Slow4 \tab 1e-21 \tab 1e-20
}
}
\section{Last updates}{


2020-11-04, DM: Added roxygen documentation
}

\examples{

# Set some arbitrary decay parameter for a dim CW-OSL measurement of quartz
name <- c("fast", "slow")
lambda <- c(2, 0.02)
n <- c(1e6, 5e7)

# Build a component table
components <- data.frame(name, lambda, n)

# Simulate the CW-OSL curve and add some signal noise
curve <- simulate_OSLcomponents(components, simulate.curve = TRUE, add.poisson.noise = TRUE)

# Perform nonlinear regression at the simulated curve
fit_results <- fit_OSLcurve(curve, K.max = 2, output.complex = TRUE)

# Plot the fitting iterations and set them into context
plot_PhotoCrosssections(fit_results)


}
\references{
Durcan, J.A., Duller, G.A.T., 2011. The fast ratio: A rapid measure for testing the dominance of the fast component in the initial OSL signal from quartz. Radiation Measurements 46, 1065–1072.

Jain, M., Murray, A.S., Bøtter-Jensen, L., 2003. Characterisation of blue-light stimulated luminescence components in different quartz samples: implications for dose measurement. Radiation Measurements 37, 441–449.

Singarayer, J.S., Bailey, R.M., 2003. Further investigations of the quartz optically stimulated luminescence components using linear modulation.
Radiation Measurements, Proceedings of the 10th international Conference on Luminescence and Electron-Spin Resonance Dating (LED 2002) 37, 451–458.
}
\seealso{
\link{fit_OSLcurve}, \link{RLum.OSL_global_fitting}
}
\author{
Dirk Mittelstraß, \email{dirk.mittelstrass@luminescence.de}

Please cite the package the following way:

Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
R package OSLdecomposition: Automated identification and separation of quartz CW-OSL signal components, \emph{in preparation}.
}
