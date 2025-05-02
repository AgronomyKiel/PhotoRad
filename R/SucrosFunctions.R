#rm(list=ls(all=TRUE)) #Löscht alle Objekte

# some functions and constants for the Sucros model ##########

# constants ################################

#conversion factor for degrees to radians
#' @export rad
rad <- pi/180.0

# number of gauss points for integration
#' @export n_gauss
n_gauss <- 5

# x-coordinates of gauss points
#' @export xgauss
xgauss <-c (0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899)

# weights of gauss points
#' @export wgauss
wgauss <-c (0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635)

# data frame for gauss points
#' @export gauss
gauss <- data.frame(xgauss, wgauss)

# transmission fraction of PAR on single leaves
#' @export transf
transf <- 0.1

# reflection fraction of PAR on single leaves
#' @export reflf
reflf <- 0.1

# scattering coefficient of leaves for visible radiation (par)
#' @export scv
scv <- transf+reflf



# simplified function for fraction of diffuse radiation based on hourly observations
# based on Spitters et al. SEPARATING THE DIFFUSE AND DIRECT COMPONENT OF GLOBAL
# RADIATION AND ITS IMPLICATIONS FOR MODELING CANOPY
# PHOTOSYNTHESIS
# PART I. COMPONENTS OF INCOMING RADIATION

#' Title
#'
#' @param Sg_S0 relation of measured to maximum solar radiation
#' @param lat latitude [°]
#' @param day day of year
#' @param hour hour of day
#'
#' @return frdif fraction of diffuse radiation
#' @export
f_frdif_hour <- function(Sg_S0, lat, day, hour){
  frdif <- pmax(ifelse(Sg_S0 < 0.22, 1,
                    ifelse(Sg_S0 >= 0.22 & (Sg_S0 < 0.35), (1- 6.4*((Sg_S0 - 0.22)^2)),
                            (1.47- 1.66*Sg_S0))),
                 0.2)

  return(frdif)
}


#' trapezoid function
#'
#' @param x independent variable
#' @param x0 lower limit where Y = fmin
#' @param x1 lower limit where Y = fmax
#' @param x2 upper limit where Y = fmax
#' @param x3 upper limit where Y = fmin
#' @param fmin lowest value of Y
#' @param fmax highest value of Y
#'
#' @returns interpolated y value of the trapezoid function
#' @export
trapez_f <- function (x, x0=0, x1=10, x2=18, x3=35, fmin = 0, fmax = 1)

# bildet folgende Funktion ab :
#
#   y |
#     |
#  fmax|             ***********
#      |           * |         | *
#      |         *   |         |   *
#      |       *     |         |     *
#      |     *       |         |       *
#  fmin ___*_________|_________|_________*_______
#          x0       x1        x2         x3
#

{
  trapez_f <- NA
if ((x >= x1) & (x <= x2)) {trapez_f <- fmax}
if ((x <= x0) | (x >= x3)) {trapez_f <- fmin}
if (x > x2) {trapez_f <- fmax - (x - x2) * (fmax - fmin) / (x3 - x2)}
if ((x > x0) & (x < x1)) {trapez_f <- fmin + (x - x0) * (fmax - fmin) / (x1 - x0)}
if (is.na(trapez_f))  {trapez_f <- 0.0}
return(trapez_f)
}


#' Title f_SucrosTransmission
#' functions for calculation of radiation transmission
#' according to the sucros model
#'
#' @param lai leaf area index
#' @param sg global radiation [J/m2/s]
#' @param lat latitude [°]
#' @param day day of year
#' @param hour hour of day
#'
#' @return transmission fraction of transmitted radiation
#' @export
f_SucrosTransmission <- function (lai, sg, lat, day, hour)
{
  dec<-f_dec(day)
  sc      <-   1370*(1+0.033*cos(2*pi*day/365));
  #  angot<-sc*dsinb;
  sinb <- f_sinbeta(lat,day,hour)
  Globmax <- sc*sinb
  parmax <- 0.5*Globmax;
  #test <- f_dec(Days)
  #plot(Days,test)
  #trans <- 0.6
  #sg <- trans*parmax
  sg_S0 <- sg/parmax
  fr_dif <- f_frdif_hour(sg_S0,lat,day,hour)
  pardf <- sg*fr_dif;
  pardr <- sg-pardf;
  sqv  <- sqrt(1.0-scv);
  vistot <- 0;

  # reflexion der horizontalen und sphärischen blattverteilung
  refh <- (1.0-sqv)/(1.0+sqv);
  refs<-refh*2.0/(1.0+1.6*sinb);

  # extinktionskoeffizienten für direkte strahlung
  kdif <- 0.8*sqv;
  kdirbl<-(0.5/sinb)*kdif/(0.8*sqv);
  kdirt<-kdirbl*sqv;
  visdf<-1:n_gauss
  vist    <- 1:n_gauss
  visd    <- 1:n_gauss
  visshd  <- 1:n_gauss
  vispp   <- 1:n_gauss
  vist    <- 1:n_gauss
  vissun  <- 1:n_gauss
  fslla   <- 1:n_gauss
  laic    <- 1:n_gauss

  # selection of depth in canopy, canopy absorption (j / m2 / s) is set to zero
  vislaic <- 0;
  vistot <- 0;
  for (layer in 1:n_gauss) {
    #layer<-1
    laic[layer]<-lai*xgauss[layer];

    # absorbed fluxes per unit leaf area: diffuse flux, total direct flux,
    # direct component of direct flux
    visdf[layer]   <-(1-refs)*pardf*kdif  *exp(-kdif  *laic[layer]);
    vist[layer]    <-(1-refs)*pardr*kdirt *exp(-kdirt *laic[layer]);
    visd[layer]    <-(1-scv) *pardr*kdirbl*exp(-kdirbl*laic[layer]);
    visshd[layer]  <-visdf[layer]+vist[layer]-visd[layer];
    vispp[layer]   <-(1-scv)*pardr/sinb;
    vissun[layer]  <- visshd[layer] + (1 - scv) * kdirbl * pardr;
    fslla[layer]   <- exp(-kdirbl*laic[layer]);
    #  fgl     <-fslla*fgrsun+(1-fslla)*fgrsh;
    #  instantaneous canopy absorption  (j / m2 / s)
    vislaic <- vislaic + (fslla[layer] * vissun[layer] + (1-fslla[layer]) * visshd[layer]) * wgauss[layer] * lai;

    #  'hourly canopy absorption (j / m2 / s)
    #  vistot  <- vistot + vislaic[layer] * wgauss [hour] * dayl  * 3600;
  }
  parabs <- vislaic;
  partrans <- sg-parabs-refs*sg
  transmission <- partrans/sg
  return(transmission)
}


