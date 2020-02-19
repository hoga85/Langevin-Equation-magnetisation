# Author: Fang-Yu Lin
# This script is used to simulation the Langevin equation.
# The Langevin equation describes the temperature dependence of magnetic magentisation. 
# It was initially proposed to described the MH curve of paramagnet.
# It has also been used to describe the MH curve of superparamagnetic materials. 
# More discussion can be found in Weaver et. al.'s article:
# J. B. Weaver, A. M. Rauwerdink, and E. W. Hansen, "Magnetic nanoparticle temperature estimation," Med. Phys., vol. 36, no. 5, pp. 1822-1829, 2009.

rm(list=ls()) 
library(stats)
library(svDialogs) # pop out box for user input
library("tcltk")

EnterPara <- T

while (EnterPara == T){
  biasedB<- as.numeric(dlgInput("What is the bias magnetic field (unit: A/m) ?",0)$res)  
  appliedB <- as.numeric(dlgInput("What is the strength of the AMF (unit: A/m) ?",2000)$res)
  appliedfreq <- as.numeric(dlgInput("What is the frequency of the AMF (unit: Hz) ?",1)$res)
  initialT <- as.numeric(dlgInput("What is the initial temperature of the MNPs (unit: K) ?",50)$res)
  finalT <- as.numeric(dlgInput("What is the final temperature of the MNPs (unit: K) ?",50)$res)
  dTemp <-  as.numeric(dlgInput("what is the incremental temperature (unit: K) ?",100)$res)
  
  Description <- paste0("Alternating Magnetic Field : ", appliedB, " A/m and ", appliedfreq," Hz", "; ",
                           "Temperature Range: between ", initialT, " (K) and",finalT," (K)", "; ",
                           "Incremental Temperature: ", dTemp, " (K)")
                     
  Answer <- tkmessageBox(message =Description,
                         icon = "question", type = "yesno", default = "yes")
  GoodPara <- (tclvalue(Answer) == "yes")
  EnterPara <-!GoodPara
}
##########################
######functions()#########
##########################
coth = function(X){
  new <- {}
  for (i in 1:length(X)){
    if (X[i] <= 710 & X[i] >= -710){
      new[i] <- cosh(X[i])/sinh(X[i])}
    else if (X[i] > 710){
      new[i] <- 1}
    else if (X[i] < -710){
      new[i] <- -1}
    else{}
  }
  return(new)}

#Domain setup
freq<-appliedfreq #Hz # signal frequency
T<- (1/freq)*2  # the duration of time, the default is two periods of the wave.
number_data = 1000
dt <- T/number_data #sampling time interval
F <-1/dt # the range of f domain
df <- 1/T # interval in f domain
t <- seq(0,T,by=dt) # the time interval in time domain
   
#CREATE OUR TIME SERIES DATA
B0 <- appliedB #A/m
B <- B0*sin(2*pi*freq*t)+biasedB
B0.all <- 10e3 # The range of the MH plot. Unit: A/m.
B.all <- B0.all*sin(2*pi*freq*t) # for MH curve.
   
#CREATE OUR FREQUENCY ARRAY
f <- 1:length(t)/T

#Parameters initialisation
i<-1
All.T <- seq(initialT,finalT,dTemp)
ratio.H3.H5 <- 0
H3<-0
H5<-0
chi.T <-0
   
#########################
# pre-assumed parameters:
#########################
   # diameter of an individual core: 130 nm 
   # (Tay, Zhi Wei et al. "In vivo tracking and quantification of 
   # inhaled aerosol using magnetic particle imaging towards inhaled 
   # therapeutic monitoring." Theranostics vol. 8,13 3676-3687. 8 Jun. 2018, 
   # doi:10.7150/thno.26608), then
   # the volume of a particle v would be:
   v <- (4/3)*((130e-9/2)^3) # unit: m3
   
   # assume the concentration of MNP suspension is 20 mg/ml
   # assume the particle density is 1.4 g/cm3
   # assume the volume of MNP suspension is 0.5 ml
   # then, the number of particle n can be estimate
   
   ###unit conversion:###
     #particle density rho:
     rho = (1.4/1000)/((1/100)^3) # unit: kg/m3
     #the weithg of all particles:
     mass_sample <- 0.5*20/1000000 # unit: kg
     #the overall volume of particles nV:
     nv <- mass_sample/rho #unit: m3
     n<- nv/v
     
   # The amplitude of Langevin Function (AL is the ratio of the magnetic energy and thermal energy):
   # To calculate the magnetic energy loss, it is important that the unit of M and H is correct so that
   # their product would have an unit of J.      
     H_tesla <- B*4*pi*1e-7 # convert the unit of A/m to T

   
   # The coil volume (coil diameter: 29.33 mm; height; 33.32 mm):
   # V_sample = (pi*(29.33e-3/2)^2)*33.32e-3 #unit: m3
   V_sample = 0.5e-6 # unit: m3
   
   # assume the Saturation magnetisation to be: 90 (A·m2/kg)
   # Convert the mass magnetisation to volume magnetisation
   M_0_mass <- 90 #unit: A·m2/kg
   M_0 <- M_0_mass*rho # unit: A/m
   
   # Boltzmann constant:
   kB <- 1.38e-23 # unit: J/K 
   
######################### 
### plot 4 figures #####   
#########################  
par(mfrow=c(2,2)) # 2x2 plots   

############################### 1st figure - MH curve ##############################
col<-0
for (Temp in All.T ){
  #index k:
  k<-Temp
  AL<- (v*M_0*H_tesla)/(4*pi*kB*Temp)    # unitless 
  M <- (nv/V_sample)*M_0*(coth(AL)-(1/AL))  # A/m  
  M <- M[2:length(M)]
  t_M <- t[2:length(t)] 
   
  H_tesla.all <- B.all*4*pi*1e-7
  AL.all<- (v*M_0*H_tesla.all)/(4*pi*kB*Temp)
  M.all <- (nv/V_sample)*M_0*(coth(AL.all)-(1/AL.all))
  M.all <- M.all[2:length(M.all)] # for MH curve 
  M.all[is.na(M.all)] <- 0
  
     #PLOTTING
  col<-col+1
  if (col ==1){
     plot(B.all[2:length(B)],M.all, type="l",lty = col,lwd=2, col = col)   
     grid()   
  
  }else{
    lines(B.all[2:length(B)],M.all, type="l",lty = col,lwd=1, col = col)   
    grid()
   
  }
     grid()
     Sys.sleep(0.1)  
     print(Temp)
}

############################### 2nd figure - M output ##############################
col<-0
for (Temp in All.T ){
  k<-Temp
  AL<- (v*M_0*H_tesla)/(4*pi*kB*Temp)    # unitless 
  M <- (nv/V_sample)*M_0*(coth(AL)-(1/AL))  # A/m  
  M <- M[2:length(M)]
  t_M <- t[2:length(t)] 
  
  H_tesla.all <- B.all*4*pi*1e-7
  AL.all<- (v*M_0*H_tesla.all)/(4*pi*kB*Temp)
  M.all <- (nv/V_sample)*M_0*(coth(AL.all)-(1/AL.all))
  M.all <- M.all[2:length(M.all)] # for MH curve 
  M.all[is.na(M.all)] <- 0
  
  #PLOTTING
  col<-col+1
  if (col ==1){
    plot(t_M,M,  type="l",lty = col,lwd=2, col = col+3,
         ylim = c(-max(M.all),max(M.all)))   
    grid()
  }else{
    lines(t_M,M,  type="l",lty = col,lwd=1, col = col,
          ylim = c(-max(M.all),max(M.all)))   
    grid()
  }
  grid()
  Sys.sleep(0.1)  
  print(Temp)
}
############################### 3rd figure - H input ##############################
col<-0
for (Temp in All.T ){
  k<-Temp
  AL<- (v*M_0*H_tesla)/(4*pi*kB*Temp)    # unitless 
  M <- (nv/V_sample)*M_0*(coth(AL)-(1/AL))  # A/m  
  M <- M[2:length(M)]
  t_M <- t[2:length(t)] 
  
  H_tesla.all <- B.all*4*pi*1e-7
  AL.all<- (v*M_0*H_tesla.all)/(4*pi*kB*Temp)
  M.all <- (nv/V_sample)*M_0*(coth(AL.all)-(1/AL.all))
  M.all <- M.all[2:length(M.all)] # for MH curve 
  M.all[is.na(M.all)] <- 0
  
  
  #PLOTTING
  col<-col+1
  if (col ==1){
    plot(B,t,  type="l",lty = col,lwd=2, col = col+1, ####
         xlim = c(-B0.all,B0.all))   
    grid()
  }else{
    lines(B,t,  type="l",lty = col,lwd=1, col = col,
          xlim = c(-B0.all,B0.all))   
    grid()
  }
  grid()
  Sys.sleep(0.1)  
  print(Temp)
  
}
############################### 4th figure - M Harmonics ##############################
col<-0
for (Temp in All.T ){
  k<-Temp
  AL<- (v*M_0*H_tesla)/(4*pi*kB*Temp)    # unitless 
  M <- (nv/V_sample)*M_0*(coth(AL)-(1/AL))  # A/m  
  M <- M[2:length(M)]
  t_M <- t[2:length(t)] 
  
  H_tesla.all <- B.all*4*pi*1e-7
  AL.all<- (v*M_0*H_tesla.all)/(4*pi*kB*Temp)
  M.all <- (nv/V_sample)*M_0*(coth(AL.all)-(1/AL.all))
  M.all <- M.all[2:length(M.all)] # for MH curve 
  M.all[is.na(M.all)] <- 0
  #FOURIER TRANSFORM WORK
  fftM <- fft(M)
  mag <- abs(fftM)*(2/length(M))

  #PLOTTING
  col<-col+1
  fft_data_length <- (length(M)/2)
  if (col ==1){

    plot(f[1:(fft_data_length)],mag[2:(fft_data_length+1)], type="l",lty = col,lwd=2, col = col,
         xlim = c(0,6*freq),ylim = c(0,2100)) ####   
  }else{

    lines(f[1:(fft_data_length)],mag[2:(fft_data_length+1)], type="l",lty = col,lwd=1, col = col,
          xlim = c(0,6*freq))   
  }
  
  grid()
  Sys.sleep(0.1)  
  print(Temp)
  
  H3[col] <- sort(mag, decreasing = TRUE)[3]
  H5[col] <- sort(mag, decreasing = TRUE)[5]
  ratio.H3.H5[col] <- sort(mag, decreasing = TRUE)[5]/ sort(mag, decreasing = TRUE)[3]   
  
}

############################### 5th figure - Ration of H3 and H5 ##############################
plot(All.T, ratio.H3.H5,type = "o")
grid()