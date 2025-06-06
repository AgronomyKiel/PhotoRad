% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PhotoFunctions.R
\name{NegExponentialr}
\alias{NegExponentialr}
\title{Negative exponential function for photosynthesis with respiration}
\usage{
NegExponentialr(Pmax, alpha, r, I)
}
\arguments{
\item{Pmax}{light saturated photosynthesis rate \link{µmol CO2 m^-2 s^-1}}

\item{alpha}{quantum yield, initial slope \link{µmol CO2 µmol^-1 photons}}

\item{r}{dark respiration rate \link{µmol CO2 m^-2 s^-1}}

\item{I}{photosynthetic photon flux density, PPFD \link{µmol photons m^-2 s^-1}}
}
\value{
single leaf net photosynthesis rate \link{µmol CO2 m^-2 s^-1}
}
\description{
Negative exponential function for photosynthesis with respiration
}
