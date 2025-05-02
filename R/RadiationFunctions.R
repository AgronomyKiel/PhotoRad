#temp <- solaR::data(helios)

# function for calculation declination of the sun (radians) as a function of Julian day
#' Title
#'
#' @param day Day of year [1-365]
#'
#' @return declination of the sun on day [radians]
#' @export
f_dec <- function(day){
  rad <- pi/180.0
  x   <- asin(sin(23.45*rad)*cos(2.*pi*(day+10)/365))
  dec <-x*-1
  return(dec)
}

# function for Daylength [h]
#' Title
#'
#' @param DAY Day of the year [1-365]
#' @param LAT Latitude [degrees]
#'
#' @return daylength [h]
#' @export
f_daylength <- function (DAY, LAT){
  RAD <- pi/180.0;

  # Zwischenwerte SINLD, COSLD und AOB }
  DEC <- f_dec(DAY)
  SINLD <-sin(RAD*LAT)*sin(DEC);
  COSLD <-cos(RAD*LAT)*cos(DEC);
  AOB   <-  SINLD/COSLD;
  # Tageslaenge (DAYL) und photoperiodische Tageslaenge (DAYLP) }
  X  <-asin(AOB);
  DAYL <- 12.0*(1+2*X/pi);
  return(DAYL)
}

# function for photoperiodic Daylength [h]
#' Title
#'
#' @param DAY day of the year [1-365]
#' @param LAT latitude [degrees]
#'
#' @return photoperiodic daylength [h]
#' @export
f_daylengthp <- function (DAY, LAT){
  RAD <- pi/180.0
  DEC <- f_dec(DAY)

  # Zwischenwerte SINLD, COSLD und AOB }

  SINLD <-sin(RAD*LAT)*sin(DEC)
  COSLD <-cos(RAD*LAT)*cos(DEC)
  AOB   <-  SINLD/COSLD;

  # Tageslaenge (DAYL) und photoperiodische Tageslaenge (DAYLP) }
  X <- asin((-sin(-4*RAD)+SINLD)/COSLD)
  DAYLP <- 12.0*(1+2*X/pi)
  return(DAYLP)
}

# function for sin beta (sun elevation)
#' Title
#'
#' @param lat Latitude [degrees]
#' @param day day of the year [1-365]
#' @param hour hour of the day [0-24]
#'
#' @return sin beta (sun elevation)
#' @export
f_sinbeta <- function(lat, day, hour){
  dec <-f_dec(day)
  rad <- pi/180.0
  sinld<-sin(rad*lat)*sin(dec)
  cosld<-cos(rad*lat)*cos(dec)
  sinbeta <- sinld+cosld*cos(2.0*pi*(hour+12.0)/24)
  return(sinbeta)
}


#' dsinbeta function
#'
#' @param LAT latitude [degrees]
#' @param DAY day of the year [1-365]
#' @param DAYL daylength [h]
#'
#' @return dsinbeta
#' @export
f_dsinbeta <- function(LAT, DAY, DAYL){
  RAD <- pi/180.0;
  DEC <- f_dec(DAY)

  # Zwischenwerte SINLD, COSLD und AOB }

  SINLD <-sin(RAD*LAT)*sin(DEC);
  COSLD <-cos(RAD*LAT)*cos(DEC);
  AOB   <-  SINLD/COSLD;
  dsinbeta <- 3600*(DAYL*SINLD+24*COSLD*sqrt(1-AOB*AOB)/pi)
  return(dsinbeta)

}


#' Fraction of leaf angle distribution
#'
#' @param chi_ratio ratio of length to width of the ellipsoid
#' @param leaf_angle angle of the leaf [degrees]
#'
#' @return fraction of leaf angle distribution
#' @export
leafAngleFractionEllipsoidal <- function(chi_ratio, leaf_angle) {
  ## Reference: Campbell, G. (1986). Extinction Coefficients for Radiation in Plant Canopies Calculated Using an Ellipsoidal Inclination Angle Distribution
  rad <- pi/180
  normalized_ellipse_area <- chi_ratio + 1.774*(chi_ratio + 1.182)^-0.733
  fraction_angle <- (2*chi_ratio^3*sin(leaf_angle*rad))/(normalized_ellipse_area*(cos(leaf_angle*rad)^2 + chi_ratio^2*sin(leaf_angle*rad)^2)^2)
  fraction_angle
}



#' calculation of chi ratio from mean leaf angle
#'
#' @param leaf_angle mean leaf angle [degrees]
#'
#' @returns chi ratio
#' @export
chiFromMeanAngle <- function(leaf_angle)	{
  ## Reference: Campbell, G. (1990). Derivation of an angle density function for canopies with ellipsoidal leaf angle distributions Agricultural and Forest Meteorology  49(3), 173-176.
  rad <- pi/180
  chi_ratio <- (9.65/(leaf_angle*rad))^(1/1.65) - 3
  chi_ratio
}


#' Function to calculate mean leaf angle from chi ratio
#'
#' @param chi ratio ratio of length to width of the ellipsoid (chi)
#'
#' @return mean leaf angle [degrees]
#' @export
MeanAngleFromChi <- function(chi)	{
  ## Reference: Campbell, G. (1990). Derivation of an angle density function for canopies with ellipsoidal leaf angle distributions Agricultural and Forest Meteorology  49(3), 173-176.
  rad <- pi/180
  angle <- 9.65*(3+chi)^-1.65
  angle <- angle*180/pi
  return(angle)
}




#' Extinction coefficient for spherical leaf angle distribution
#'
#' @param sun_angle angle of the sun [degrees]
#'
#' @returns extinction coefficient for spherical leaf angle distribution
#' @export
k_spherical <- function(sun_angle) {
  theta <- sun_angle * (pi / 180)
  # Spherical LAD: k = 0.5 / cos(theta)
  0.5 / cos(theta)
}



#' Extinction coefficient for cylindrical leaf angle distribution
#'
#' @param sun_angle angle of the sun [degrees]
#'
#' @returns extinction coefficient for cylindrical leaf angle distribution
#' @export
k_cylindrical <- function(sun_angle) {
  theta <- sun_angle * (pi / 180)
  # Cylindrical LAD: k = 2 / (pi * tan(theta))
  2 / (pi * tan(theta))
}



#' Calculation of extinction coefficient for ellipsoidal leaf angle distribution
#'
#' @param chi_ratio ratio of length to width of the ellipsoid (chi)
#' @param sun_angle angle of the sun [degrees]
#'
#' @return extinction coefficient for ellipsoidal leaf angle distribution
#' @export
extinctionCoefficientEllipsoidal <- function(chi_ratio, sun_angle) {
  ## Reference: Campbell, G. S., & Norman, J. (2012). An introduction to environmental biophysics. Springer Science & Business Media.
  rad <- pi/180
  sun_angle <- 90-sun_angle
  normalized_ellipse_area <- chi_ratio + 1.774*(chi_ratio + 1.182)^-0.733
  extinction <- sqrt(chi_ratio^2 + tan(sun_angle*rad)^2)/normalized_ellipse_area
  extinction
}



#' Calculation of extinction coefficient for diffuse radiation
#' with ellipsoidal leaf angle distribution
#'
#' @param lai leaf area index [dimensionless]
#' @param chi_ratio ratio of length to width of the ellipsoid (chi)
#'
#' @return extinction coefficient for diffuse radiation
#' @export
KdifEllipsoidal <- function(lai, chi_ratio){
  A<- 0.65
  B<- 1.9
  kd_inf <- 2/pi*atan(chi_ratio)
  LAI_exp_A <- lai^A
  kd <- (kd_inf*LAI_exp_A+B)/(LAI_exp_A+B)
  return(kd)
}





#' Title sun elevation angle
#'
#' @param solar.hour Solar hour [0-24]
#' @param latitude Latitude [degrees]
#' @param day.of.year Day of the year [1-365]
#'
#' @return sun elevation angle [degrees]
#' @export
Elevation <- function(solar.hour, latitude, day.of.year) {
  rad <- (pi/180)
  declination <- asin(sin(-23.45*rad)*cos((360*(day.of.year + 10)/365)*rad))/rad
  a <- sin(latitude*rad)*sin(declination*rad)
  b <- cos(latitude*rad)*cos(declination*rad)
  sinus.elevation <- a + b*cos(15*((solar.hour - 12))*rad)
  elevation <- asin(sinus.elevation)/rad
  return(elevation)
}


# simplified function for fraction of diffuse radiation based on hourly observations
# based on Spitters et al. SEPARATING THE DIFFUSE AND DIRECT COMPONENT OF GLOBAL
# RADIATION AND ITS IMPLICATIONS FOR MODELING CANOPY
# PHOTOSYNTHESIS
# PART I. COMPONENTS OF INCOMING RADIATION

#' Estimation of fraction of diffuse radiation on hourly basis
#'
#' @param Sg_S0 ratio of actual to maximum global radiation [-]
#'
#' @return fraction of diffuse radiation [-]
#' @export
f_frdif_hour <- function(Sg_S0){
  frdif <- pmax(ifelse(Sg_S0 < 0.22, 1,
                       ifelse(Sg_S0 >= 0.22 & (Sg_S0 < 0.35), (1- 6.4*((Sg_S0 - 0.22)^2)),
                              (1.47- 1.66*Sg_S0))),
                0.2)

  return(frdif)
}


#' fraction of diffuse radiation based on daily observations
#'
#' @param ATMTR atmospheric transmission
#'
#' @return fraction of diffuse radiation [-]
#' @export
f_frdif_day <- function(ATMTR){

  frdif <- ifelse(ATMTR > 0.75, 0.23,
                  ifelse((ATMTR<=0.75)&(ATMTR>0.35),1.33-1.46*ATMTR,
                         ifelse((ATMTR <= 0.35) & (ATMTR > 0.07), 1-2.3*(ATMTR-0.07)*(ATMTR-0.07),1)))
  return(frdif)
}



#' calculation of radiation intensity on a sloped plane from sunAngle
#' and Radiation intensity of a flat surface
#'
#' @param BeamFlat radiation intensity on a flat surface [W/m²]
#' @param SunAngle angle of the sun [degrees]
#' @param LeafAngle angle of the leaf [degrees]
#'
#' @return radiation intensity on a sloped leaf plane [W/m²]
#' @export
BeamOnLeaf <- function (BeamFlat, SunAngle, LeafAngle){
# calculation of radiation intensity on a sloped plane from sunAngle and Radiation intensity of
  # a flat surface
  rad <- pi/180
  Beam <- BeamFlat*sin(SunAngle*rad+LeafAngle*rad)/sin(SunAngle*rad)

}


#' calculation of radiation intensity on a sloped plane from sunAngle and Radiation intensity of
#' a flat surface
#'
#' @param BeamFlat radiation intensity on a flat surface [W/m²]
#' @param SunAngle angle of the sun [degrees]
#' @param LeafAngle angle of the leaf [degrees]
#' @param SunAzimuthAngle Azimuth angle of the sun [degrees]
#' @param LeafAzimuthAngle Azimuth angle of the leaf [degrees]
#'
#' @return radiation intensity on a sloped leaf plane [W/m²]
#' @export
BeamOnLeafAzimuth <- function (BeamFlat, SunAngle, LeafAngle, SunAzimuthAngle, LeafAzimuthAngle){
  #
  rad <- pi/180
  SunAngle <- SunAngle*rad
  LeafAngle <- LeafAngle*rad
  SunAzimuthAngle <- SunAzimuthAngle*rad
  LeafAzimuthAngle <- LeafAngle*rad
  Beam <- BeamFlat*(cos(SunAngle)*sin(LeafAngle)*cos(LeafAzimuthAngle-SunAzimuthAngle)+
                        sin(SunAngle)*cos(LeafAngle))

}

#' Estimation of net radiation
#'
#' @param albedo albedo of the surface [-]
#' @param R_s global radiation [W/m²]
#' @param R_sclear clear sky radiation [W/m²]
#' @param Airtemp air temperature [°C]
#' @param e_a actual vapor pressure [hPa]
#'
#' @return net radiation [W/m²]
#' @export
r_net <- function(albedo=0.77, R_s, R_sclear, Airtemp=20, e_a)
{
  sigma <- 5.67e-8
  abs0 <- 273.15
  emissivity <- 0.98
  T_a <- Airtemp+abs0
  R_n <- albedo*R_s+(R_s/R_sclear)*(1.24*((10*e_a)/T_a)^(1/7)-1)*emissivity*sigma*T_a^4
  r_net <- R_n
}

