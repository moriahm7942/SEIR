##################################################################################
##################################################################################
# an R script to compare the predictions of a simple SIR model to influenza epidemic
# data
#
# Author:  Sherry Towers
#          smtowers@asu.edu
# Created: Dec 6, 2012
#
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
##################################################################################

##################################################################################
# par(pch=20) sets a solid round dot style for the plots
# The chron package contains utilities to calculate dates
##################################################################################
rm(list = ls(all = TRUE))  # resets R to fresh
require("chron")
require("sfsmisc")
par(pch=20)  
source("sir_func.R")

##################################################################################
# read in 2007 influenza surveillance data (confirmed cases) for
# the Midwest (CDC region 5) from
# http://www.cdc.gov/flu/weekly/regions2007-2008/datafinalHHS/whoregX.htm
# where X=5 is midwest
# X=1 is northeast
# X=2 is NY and NJ
# X=3 are eastern seabord states like PA, DE, etc
#
# the weeks are number of weeks relative to Jan 1, 2007
# week 1 in 2007 ended Jan 6, 2007
##################################################################################
adat = read.table("midwest_influenza_2007_to_2008.dat",header=T,sep=",")
cat("\n")
cat("The data file contains: ",names(adat),"\n")
cat("\n")

##################################################################################
# The CDC data is weekly, with the date of each point corresponding to the
# date of the end of the week over which the data were collected.
# Let's convert these dates to time in days, relative to Jan 1, 2007
# We will be using this vector of dates, vtime_data, to obtain the model estimates
# of the incidence at that time.
# adat$week is relative to Jan 1, 2007, with week #1 occuring the first week in
# January, 2007.
##################################################################################
adat$time_in_days_rel_jan_1_2007 = julian(1,6,2007)+(adat$week-1)*7-julian(1,1,2007)

##################################################################################
# Specifically for this data:
# make sure we are far enough into the season that there is at least one case per week
##################################################################################
adat=subset(adat,week>=47)  
incidence_observed = adat$B
times_of_observed = adat$time_in_days_rel_jan_1_2007
time_binning = min(diff(times_of_observed))


##################################################################################
##################################################################################
# set up the model parameters
# npop is approximately the population of IL IN MI MN OH WI (CDC region 5)
# I_0 is the initial number infected
# R_0 is the inital number recovered and immune
# S_0 is the susceptibles
#
# 1/gamma is the average recovery period of the disease
# R0      is the reproduction number
# t0      is the time-of-introduction of the disease to the population, measured
#         in days from Jan 1, 2007
#
# For the SIR model, R0=beta/gamma, thus given our hypotheses for gamma and R0,
# we calculate beta as beta=R0*gamma (note that beta and gamma appear in the model
# equations, not gamma and R0, which is why we need to calculate beta given R0
# and gamma).
##################################################################################
npop = 52000000 
I_0 = 1      
R_0 = 0      
S_0 = npop-I_0-R_0

gamma = 1/3  
#R0 = 1.24    
#t0 = 233     
R0 = 1.35    
t0 = 310     
beta  = R0*gamma

vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)
   
##################################################################################
# We get the model solution for all days from t0 to the last week of the
# data time series.  If t0 is greater than the minimum date in the data time
# series, we need to print out a warning, because the time of introduction had
# to be before we actually started observing cases in the data!
##################################################################################
tmin = t0
tmax = max(times_of_observed)

if (tmin>(min(times_of_observed)-time_binning)){
  cat("\n")
  cat("**************************************************************************\n")
  cat("**************************************************************************\n")
  cat("The time-of-introduction is _after_ the first cases appeared!",t0,min(times_of_observed)-time_binning,"\n")
  cat("**************************************************************************\n")
  cat("**************************************************************************\n")
  cat("\n")
}
tmin = min(t0,min(times_of_observed)-time_binning)

vt = seq(tmin,tmax)

##################################################################################
# Now solve the system of differential equations numerically with lsoda in the 
# deSolve package.  Put the results in solved_model
# The derivative_calc_func for the SIR model is in the sir_func.R script
##################################################################################
solved_model = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters)) 

##################################################################################
# the B influenza data is incidence, not prevalence (actually, the B influenza
# data is the true incidence times the fraction that the CDC actually confirms)
#
# The incidence over some time step is the difference in S over that time
# step (in a model with no births or deaths, immigration or emigration, or 
# recurring susceptibility).
#
# The time step are derived from the dates at the end of the week of each
# data point (vtime_data)
#
# solved_model$time%in%vtime_data returns the indices of the elements of simodel$time that
# are also found in vtime_data
##################################################################################
tmin_data = min(times_of_observed)-time_binning
tmax_data = max(times_of_observed)
vtime_data = seq(tmin_data,tmax_data,time_binning)

susceptible_predicted = solved_model$S[solved_model$time%in%vtime_data]  
incidence_predicted = -diff(susceptible_predicted)

##################################################################################
# from the model estimate of the incidence and the data, we can 
# estimate the fraction of cases that were confirmed
##################################################################################
frac_confirmed = sum(incidence_observed)/sum(incidence_predicted)
cat("\n")
cat("The estimated fraction of confirmed cases among all infections is = ",frac_confirmed,"\n")
cat("\n")

##################################################################################
# normalize the model prediction so area under curve
# equals the sum of the data incidence
##################################################################################
incidence_predicted = incidence_predicted*frac_confirmed 
   
##################################################################################
# now let's overlay the model on the data, and calculate the least-squares
# statistic that compares the data to this model calculated
# under a particular hypothesis of R0 and t0
##################################################################################
if (length(incidence_predicted)>1
   &!is.na(sum(incidence_predicted))){

   if (min(incidence_predicted)>0
      &length(incidence_predicted)==length(incidence_observed)){

      ######################################################################## 
      # calculate the least squares statistic
      ######################################################################## 
      least_squares = sum((incidence_predicted-incidence_observed)^2)

      ######################################################################## 
      # plot the results
      # cex is the point size
      ######################################################################## 
      mult.fig(1,main="Confirmed B influenza cases, Midwest region, 2007-2008 season")
      ymax = max(c(incidence_observed,incidence_predicted))
      plot(times_of_observed
          ,incidence_observed
          ,ylim=c(0,1.2*ymax)
          ,xlab="Time, in weeks relative to Jan 1, 2007"
          ,ylab="Incidence"
          ,cex=2)
      lines(times_of_observed,incidence_predicted,col=2,lwd=5) # overlay the model

      ######################################################################## 
      # show the distance between the model and data for each point
      ######################################################################## 
      for (iind in 1:length(times_of_observed)){
         arrows(times_of_observed[iind]
               ,incidence_predicted[iind]
               ,times_of_observed[iind]
               ,incidence_observed[iind]
               ,code=3
               ,lwd=1
               ,length=0.10
               ,col=4
               ,angle=20)
      }
      legend("topleft"
            ,legend=c("Data","SIR model prediction","Point-by-point distance between data and model")
            ,col=c(1,2,4)
            ,lwd=3
            ,bty="n")

      ######################################################################## 
      # print the least squares info on the plot
      # first overlay a blank plot with axis limits 0 to 1 on both x and y
      # this way we can always have the text appear at the same relative place on the plot
      ######################################################################## 
      par(new=T)
      plot(c(0,1),c(0,1),axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",col=0)

      xtext = 0.03
      cex_text = 0.75
      text(xtext,0.65,paste("Model initial conditions and parameters:"),adj=0,cex=cex_text)
      text(xtext,0.6,paste(" population = ",npop,sep=""),adj=0,cex=cex_text)
      text(xtext,0.55,paste(" I0 = ",I_0," and rest of population susceptible",sep=""),adj=0,cex=cex_text)
      text(xtext,0.5,paste(" 1/gamma  = ",1/gamma," days",sep=""),adj=0,cex=cex_text)
      text(xtext,0.45,paste(" t0  = week ",round(t0/7),sep=""),adj=0,cex=cex_text)
      text(xtext,0.4,paste(" R0 = ",R0,sep=""),adj=0,cex=cex_text)
      text(xtext,0.30,paste(" The least-squares \n statistic is = ",round(least_squares,1),sep=""),adj=0,cex=cex_text)

   } # end check that the length of the incidence vector matches the length of the data vector
} # end check that the incidence vector actually has some data in it




