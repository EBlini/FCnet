% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFC.R
\name{plotFC}
\alias{plotFC}
\title{plot a Functional Connectivity matrix}
\usage{
plotFC(
  FCmatrix,
  style = c("lower.tri", "full"),
  col = optionsFCnet("colorPaletteDefault"),
  limit = NULL,
  network_definition = NULL,
  plot_labels = F,
  colorNuances = optionsFCnet("colorNuances")
)
}
\arguments{
\item{FCmatrix}{The FC matrix to plot}

\item{style}{Whether to plot only the lower triangle (and diagonal) or the entire FC matrix.}

\item{col}{The color palette used to plot values (i.e. a vector of colors). The vector is often divergent, es. \code{c("red", "white", "blue")}.
Defaults to \code{optionsFCnet("colorPaletteDefault")}.}

\item{limit}{The limits for the FC values. Defaults to \code{NULL} and automatically adapts to the data range.}

\item{network_definition}{A character vector specifying to which FC network the ROI belongs to. If provided, draws vertical and horizontal lines visually separating the networks.}

\item{plot_labels}{if \code{TRUE} annotates the network names defined in \code{network_definition}. Pretty results not warranted.}

\item{colorNuances}{Number of nuances along the provided colorPalette.}
}
\value{
A \code{ggplot2} object which can be customized further.
}
\description{
This function reads a (square) matrix and plots it using \code{ggplot2}.
}
