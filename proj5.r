## Code & comments produced by Evangelos Antypas :: s2449453
##-----------------------------------------------------------------------------
## Summary :: To be completed
##-----------------------------------------------------------------------------

## Fixing the working directory
## This line will be commented out in the final version
library(rjags)
setwd("C:\\Users\\Vaggelis Antypas\\SP_Practical5")

##-----------------------------------------------------------------------------

## Loading the data 

## death1722 contains :: 
## - deaths ::  the number of deaths that week
## - week ::  the week since the start of 2017, where 2020 starts in week 157
## - d :: the mortality rate modifier, dj , value for that week

death1722 <- read.table("death1722uk.dat", header = TRUE)

## Extracting the columns we are going to use by name

deaths <- death1722[, "deaths"]
weeks  <- death1722[, "week"]
d <- death1722[, "d"]

##-----------------------------------------------------------------------------

## lt1720 contains ::
## - age :: the age classes 0:0-1 years, 1: 1-2 years, etc.
## - fpop17 and mpop17 :: the female and male populations in each 1 
## year age class at the start of 2017
## - fpop20 and mpop20 :: the same for the start of 2020
## - mf and mm :: female and male annual death rates for each 1-year age 
## band, computed from 2017-19 data

lt1722 <- read.table("lt1720uk.dat", header = TRUE)

## Extracting the columns we are going to use by name

ages <- lt1722[, "age"]
fpop17 <- lt1722[, "fpop17"]
mpop17 <- lt1722[, "mpop17"]
fpop20 <- lt1722[, "fpop20"]
mpop20 <- lt1722[, "mpop20"]
mf <- lt1722[, "mf"]
mm <- lt1722[,"mm"]

##-----------------------------------------------------------------------------

## Defining the death_pred function

death_pred <- function(Nm,Nf,mm,mf,d){
  
  ## Inputs :: 
  ## - Nm,Nf :: vectors containing the male and female populations for each age 
  ## class (0-100) respectively
  ## - mm,fm :: the instantaneous per capita death rate per year for males and 
  ## females respectively
  ## - d :: adjusts for seasonality (per week) in the death rates
  ## (irrespective of gender)
  ##---------------------------------------------------------------------------
  
  ## Outputs :: 
  ## - total_deaths :: a vector containing the total number of  deaths for 
  ## every week supplied by the user 
  ##---------------------------------------------------------------------------
  
  ## Initializations
  
  weeks <- length(d)
  ages <- length(Nm)
  deaths <- rep(0,weeks)
  
  ## Male population
  
  Dm <- rep(0,ages)
  qm <- rep(0,ages)
  Nm_star <- rep(0, ages)
  Nm_plus <- rep(0, ages)
  
  N0m_star <- Nm[1] 
  
  ## Female population
  
  Df <- rep(0,ages)
  qf <- rep(0,ages)
  Nf_star <- rep(0, ages)
  Nf_plus <- rep(0, ages)
  
  N0f_star <- Nf[1] 
  
  for(j in c(1:weeks)){
    
    
    for(i in c(1:ages)){
      
      ## Male and Female per capita mortality rates
      
      qm[i] <- 1 - exp(-mm[i]/52)
      qf[i] <- 1 - exp(-mf[i]/52)
      
      ## Male and Female Weekly deaths
      
      Dm[i] <- 0.9885*d[j]*qm[i]*Nm[i]
      Df[i] <- 0.9885*d[j]*qf[i]*Nf[i]
      
      ## We update the populations by taking out the people who died, for males
      ## and females respectively 
      
      Nm_star[i] <- Nm[i] - Dm[i]
      Nf_star[i] <- Nf[i] - Df[i]
      
      ## This procedure applies aging by a week to each age class
      
      if(i==1){ ## For the first age class we have initialized with N[1] which 
                ## is assumed constant 
        
        Nm_plus[i] <- (51/52)*Nm_star[i] + (1/52)*N0m_star
        Nf_plus[i] <- (51/52)*Nf_star[i] + (1/52)*N0f_star
        
      } else{
        
        Nm_plus[i] <- (51/52)*Nm_star[i] + (1/52)*Nm_star[i-1]
        Nf_plus[i] <- (51/52)*Nf_star[i] + (1/52)*Nf_star[i-1]
      }
      
      
    }
    
    ## We update the populations after weekly deaths and ageing 
    
    Nm <- Nm_plus
    Nf <- Nf_plus
    
    ## we sum male and female deaths to obtain to total weekly deaths
    
    deaths[j] <- sum(Dm + Df)
  }
  
  return(deaths)
  
} 

## Excess deaths from the start of the year 2020 to the end of the data 
## 2020 starts from week 157 

preds <- death_pred(mpop20,fpop20,mm,mf,d[157:305])
excess_deaths <- sum(deaths[157:305]) - sum(preds)


## Excess deaths just for the year 2020, a year consists of 48 weeks so the 
## year will end at week 157 + 48 = 205 

preds_2020 <- death_pred(mpop20,fpop20,mm,mf,d[157:208])
excess_deaths_2020 <-  sum(deaths[157:208]) - sum(preds_2020)

## Plotting observed deaths from the start of 2020 vs weeks, along with the 
## excess deaths for the year 2020

plot(weeks[157:305], deaths[157:305],
     main = paste(" Actual Deaths from 2020 vs Weeks, 
                  Total Excess deaths from 2020::",floor(excess_deaths),
                  "//Excess deaths in 2020::",floor(excess_deaths_2020)),
     xlab = "Weeks", ylab = "Deaths", pch = 19,
     ylim = c(0,max(deaths[157:305])))

## Adding excess deaths on the graph 

lines(weeks[157:305], preds, type = "l", col = "red", lwd=2.0)

## Adding a legend 
  
legend("topleft", legend=c("Actual Deaths", "Predicted Deaths"),
       col=c("black","red"), lty = c(NA, 1), 
       pch = c(19, NA), cex=0.8, bg = "lightgreen")
  
tot_preds <- death_pred(mpop20,fpop20,mm,mf,d[157:305])
vec_excess_deaths <- deaths[157:305] - tot_preds

cumsum_excess_deaths <- cumsum(vec_excess_deaths)

plot(weeks[157:305],cumsum_excess_deaths,main = "Cummulative excess deaths vs Weeks", 
     xlab = "Weeks", ylab = "Excess deaths", pch = 19)
  
  
legend("topleft", legend=c("Cummulative excess deaths"),
       col=c("black"), 
       pch = c(19), cex=0.8, bg = "lightgreen")
  
  
  
## JAGS model

## library(rjags)
model <- jags.model("model.jags",data=list(x=deaths,N=length(deaths)))
  
  
  
  
  
  
  
  
  
  
  