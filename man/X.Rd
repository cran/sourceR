% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{X}
\alias{X}
\title{Constructs source data}
\usage{
X(data, x, type, source, time = NULL)
}
\arguments{
\item{data}{long-format data.frame containing source data}

\item{x}{character string giving name of source counts column in data}

\item{type}{character string giving name of type column in data}

\item{source}{character string giving name of source column in data}

\item{time}{optional column denoting times of source observation}
}
\value{
A X source data structure for use in sourceR models
}
\description{
The X constructor function returns an R6 X class
which feeds source data into sourceR models.
}
