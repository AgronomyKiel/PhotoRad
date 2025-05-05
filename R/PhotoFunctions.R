
#' Blackman function for photosynthesis
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
Blackman <- function(Pmax, alpha, I){
  result <- ifelse(I*alpha>Pmax,Pmax,alpha*I)
  return(result)
}

#' Blackman function for photosynthisis with respiration
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param r dark respiration rate [µmol CO2 m^-2 s^-1]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf net photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
Blackmanr <- function(Pmax, alpha, r, I){
  result <- ifelse(I*alpha>Pmax,Pmax,alpha*I)-r
  return(result)
}


#' Rectangular hyperbola function for photosynthesis
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
NegExponential<- function(Pmax, alpha, I){
  result <- Pmax*(1-exp(-alpha*I))
  return(result)
}

#' Negative exponential function for photosynthesis with respiration
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param r dark respiration rate [µmol CO2 m^-2 s^-1]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf net photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
NegExponentialr<- function(Pmax, alpha,r, I){
  result <- Pmax*(1-exp(-alpha*I))-r
  return(result)
}



#' Rectangular hyperbola function for photosynthesis
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
RectHyp <- function(Pmax, alpha, I){
  result <-alpha*Pmax*I/((alpha*I)+Pmax)
  return(result)
}


#' Rectangular hyperbola function for photosynthesis with respiration
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param r dark respiration rate [µmol CO2 m^-2 s^-1]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf net photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
RectHypR <- function(Pmax, alpha, r, I){
  result <- alpha*Pmax*I/((alpha*I)+Pmax)-r
  return(result)
}



#' Non-rectangular hyperbola function for photosynthesis without respiration
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param theta curvature parameter [dimensionless]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf net photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
NonRectHyp <- function(Pmax, alpha, theta, I){
  result <- (alpha*I+Pmax-sqrt((alpha*I+Pmax)^2-4*theta*alpha*I*Pmax))/(2*theta)
  return(result)
}




#' Non-rectangular hyperbola function for photosynthesis with respiration
#'
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param r dark respiration rate [µmol CO2 m^-2 s^-1]
#' @param theta curvature parameter [dimensionless]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns single leaf net photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
NonRectHypR <- function(Pmax, alpha, r, theta, I){
  result <- (alpha*I+Pmax-sqrt((alpha*I+Pmax)^2-4*theta*alpha*I*Pmax))/(2*theta)-r
  return(result)
}



#' Canopy Integral of rectangular hyperbola function for photosynthesis with constant Pmax
#'
#' @param LAI leaf area index [dimensionless]
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param k extinction coefficient [dimensionless]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns canopy-integrated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
RectHypInt <- function(LAI, Pmax, alpha,  k, I) {
  result <- Pmax/k*log((alpha*k*I+Pmax)/(alpha*k*I*exp(-k*LAI)+Pmax))
  return(result)
}


#' Canopy Integral of rectangular hyperbola function for photosynthesis with constant Pmax and respiration
#'
#' @param LAI leaf area index [dimensionless]
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param r dark respiration rate per unit leaf area [µmol CO2 m^-2 s^-1]
#' @param k extinction coefficient [dimensionless]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns canopy-integrated net photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
RectHypIntR <- function(LAI, Pmax, alpha, r, k, I) {
  result <- Pmax/k*log((alpha*k*I+Pmax)/(alpha*k*I*exp(-k*LAI)+Pmax))-r*LAI
  return(result)
}


#' Canopy integral of rectangular hyperbola function for photosynthesis with decrease Pmax of Pmax
#'  and respiration
#'
#' @param LAI leaf area index [dimensionless]
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param k extinction coefficient [dimensionless]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns canopy-integrated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
RectHypIntPmaxdec <- function(LAI, Pmax, alpha, k, I) {
  result <- alpha*I*Pmax*(1-exp(-k*LAI))/(alpha*k*I+Pmax)
  return(result)
}




#' Canopy Integral of rectangular hyperbola function for photosynthesis with decrease Pmax of Pmax
#' and respiration
#'
#' @param LAI leaf area index [dimensionless]
#' @param Pmax light saturated photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @param alpha quantum yield, initial slope [µmol CO2 µmol^-1 photons]
#' @param r dark respiration rate per unit leaf area [µmol CO2 m^-2 s^-1]
#' @param k extinction coefficient [dimensionless]
#' @param I photosynthetic photon flux density, PPFD [µmol photons m^-2 s^-1]
#'
#' @returns canopy-integrated net photosynthesis rate [µmol CO2 m^-2 s^-1]
#' @export
RectHypIntPmaxdecR <- function(LAI, Pmax, alpha, r, k, I) {
  result <- alpha*I*Pmax*(1-exp(-k*LAI))/(alpha*k*I+Pmax)-(-r/k*exp(-k*LAI)+r/k)
  return(result)
}



#' Utility function for transfer of a daily average to a 16 hour sinusoidal function
#'
#' @param hour hour of the day [0-24]
#'
#' @returns the actual value of the sinusoidal function for that hour
#' @export
sinusf <- function (hour=12)

{
  output <-  pmax(0,1.64221194*(0.5+sin(pi*((hour+18)/12))))
  #                   1.64221194*(0.5+sin(pi*((Stunden+18)/12)))
  return(output)
}



