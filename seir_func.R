##################################################################################
# An R script to solve ODE's of an SEIR model 
# Modified from:
# http://www.sherrytowers.com/sir_func.R
#
# Author: Sherry Towers
#         smtowers@asu.edu
#
# Created: Dec 1st, 2012
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors.
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
#
# In order to use this script in R, you must first have the deSolve library
# installed.  To do this, type in the R command line:
#    install.packages("deSolve")
# then pick a mirror site located close to you.
##################################################################################

require("deSolve")

##################################################################################
##################################################################################
# This is a function which, given a value of S, E,I, and R at time t
# calculates the time derivatives of S E I and R
#
# The function gets called over and over again by functions in the R deSolve
# library to do the numerical integration over many time steps to solve
# the system of ODE's
#
# t is the current time at a particular step
# The vector x contains the current values in each compartment
# The list object vparameters contains the parameters of the model, like the
# recovery period, gamma, and the transmission rate, beta (in this case)
# In the main program, the vparameters list object is filled with the parameter
# values and their names.  
#
# This function gets passed to the functions in the deSolve package
#
##################################################################################
derivative_calc_func=function(t, x, vparameters){
   ###############################################################################
   # The vector x is the same length as the number of compartments.  In your main
   # program you identify which compartment corresponds to which element of x.
   # You need to make sure that these are in the same order. 
   ###############################################################################
   S = x[1]  
   E = x[2]  
   I = x[3]
   R = x[4]  

   ###############################################################################
   # vparameters is a list object, filled in the main program, and passed
   # Keep this next line the same when you are writing your own function to
   # solve a system of ODE's
   ###############################################################################
   with(as.list(vparameters),{

      ############################################################################
      # calculate the population size, which for a simple SIR model is
      # npop = S+I+R
      # we will need this to calculate our SIR model derivatives below
      ############################################################################
      npop = S+E+I+R   

      ############################################################################
      # Now give the expressions for the derivatives of S, I, and R wrt time
      # these come from the model equations.  The following equations are for
      # an SEIR model.  When you write your own function, replace these with
      # your model equations
      ############################################################################
      dS = -beta*S*I/npop     
      dE = beta*S*I/npop-sigma*E    
      dI = sigma*E-gamma*I  
      dR = gamma*I              

      ############################################################################
      # vout is an output vector that contains the values of the derivates
      # the compartments in the vector on the RHS need to be in the same order as the
      # compartments used to fill the x vector!
      ############################################################################
      vout = c(dS,dE,dI,dR)
      list(vout)
   })
}

##################################################################################
##################################################################################
##################################################################################
# this is the same as the above function, except now it includes births and
# deaths (both with rate mu) in the model equations
##################################################################################
derivative_calc_func_with_demographics=function(t, x, vparameters){
   S = x[1]  
   E = x[2]  
   I = x[3]
   R = x[4] 

   with(as.list(vparameters),{
      npop = S+E+I+R   
      dS = -delta*npop+beta*S*I+lambda*S     
      dE = -beta*S*I+gamma*E+lambda*E     
      dI = -gamma*E+mu*I+lambda*I  
      dR = -mu*I+lambda*R                      
      out = c(dS,dE,dI,dR)
      list(out)
   })
}

##################################################################################
##################################################################################
# this is the same as the derivative_calc_function, except now this involves
# calculating the derivatives of an SIR model with vaccination
# Rvac is the vaccinated (and now immune) compartment
# The vaccination begins at time_vaccination_begins, and ends at 
# time_vaccination_ends
##################################################################################
derivative_calc_func_with_vaccination=function(t, x, vparameters){
   S = x[1]  
   E = x[2]  
   I = x[3]
   R = x[4]
   Rvac = x[5]  

   with(as.list(vparameters),{
      npop = S+E+I+R+Rvac   
      dS    = mu*N-mu*S-beta*S*I/npop-nu*S            
      dE    = +beta*S*I/npop - mu*E-sigma*E
      dI    = +sigma*E-mu*I-gamma*I  
      dR    = +gamma*I-mu*R+nu*S                  
      dRvac = 0
      if (t>=time_vaccination_begins&t<=time_vaccination_ends){
         dS    = dS - rho*S
         dRvac = +rho*S
      }
      out = c(dS,dE,dI,dR,dRvac)
      list(out)
   })
}



