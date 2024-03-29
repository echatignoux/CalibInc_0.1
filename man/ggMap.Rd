% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggMap.R
\name{ggMap}
\alias{ggMap}
\title{Plot variables values on a map}
\usage{
ggMap(var, data = NULL, map = NULL, limits = NULL, breaks = NULL,
  legend = list(breaks = NULL, density = F, title = NULL, pretty = T,
  placement = NULL), color = "Greys", rev = FALSE, na.action = na.pass,
  ns = NULL, path = TRUE)
}
\arguments{
\item{var}{Variable name from the data \code{data} to map, or
vector sorted in the order of the map spatial units}

\item{data}{Data with column to plot. Must contain an \code{id} variable
that identifies the district (jointure variable with \code{map} data)}

\item{map}{Map shape ("polygg" type, set to ShapDep if NULL)}

\item{limits}{Set limits if a continuous map is chosen}

\item{breaks}{Values to cut \code{var} into classes}

\item{legend}{List of legend parameters \describe{
\item{breaks}{Breaks in the legend}
\item{density}{Set to \code{TRUE} to represent the legend as the density distribution
of the observed values}
\item{title}{Legend title (variable name otherwise)}
\item{pretty}{If \code{var} is cut in classe with \code{cut} function, convert classes names to pretty labels}
\item{placement}{Lengend placement (x0,x1,y0,y1) coordinates,
in proportion of the plot zone} }}

\item{color}{Name of the brewer palette
\code{\link[RColorBrewer]{brewer.pal}}}

\item{rev}{Reverse color palette order (\code{TRUE,FALSE})}

\item{na.action}{What to do with missing values}

\item{ns}{Variable name in the data set that identifies non significant areas (\code{TRUE} when
non significant). In that case, the corresponding geographical areas are hatched.}

\item{path}{Draw geographical boundaries (default to \code{TRUE})}
}
\value{
A ggplot2 graph
}
\description{
\code{ggMap} create a map from a data frame with observations on spatial units.
It allows to represent values in continuous (with possible lower and upper bounds)
or categorial scale, and to hatch non significant areas.
}
\examples{
## Represent crude hospitalisation rate in France by district for LP cancer in men, 2007-2011
library(tidyverse)
data(lopm.alldist)
crude.rate<-lopm.alldist\%>\%group_by(dist)\%>\%summarise(rate=100000*sum(H)/sum(py))\%>\%mutate(id=dist)
ggMap(var=rate,data=crude.rate,color= "YlOrRd",legend=list(title="Hosp. rate\\nfor 100 000"))+
ggtitle("Crude hospitalisation rate for LOP, men, 2007-2011")
}
\author{
Edouard Chatignoux
}
