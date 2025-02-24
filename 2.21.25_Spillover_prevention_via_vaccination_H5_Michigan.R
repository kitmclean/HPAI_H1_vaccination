

### Basic H5 spillover ODE for Huron County, MI (2024) 

library(deSolve)
library(tidyr)
library(ggplot2)
library(dplyr)

#setwd("C:/Users/kathe/Dropbox (University of Michigan)/Manuscripts/HPAI in cattle")


# Estimating dairy worker population in a county requires us to correct for undocumented workers. However, Michigan does not track numbers of dairy workers,
# as these are not seasonal positions. 

# Approximately 44% of undocumented agricultural workers do not self-identify (https://www.farmworkerjustice.org/wp-content/uploads/2022/06/NAWS-data-fact-sheet-FINAL.docx-3.pdf)

# Lewis, B., Martinez, R., & Coronado, J. (2017). Farmworkers in Michigan. JSRI Report Report No. 29.
# East Lansing, MI: The Julian Samora Research Institute, Michigan State University. 

# PopSize = 1217 * .56 * 5 # = (No. of farms in Huron County * % that are livestock-based * national avg # of employees per dairy) 

# Data sources: 
#No. farms and % livestock-based data from: https://www.nass.usda.gov/Publications/AgCensus/2022/Online_Resources/County_Profiles/Michigan/cp26063.pdf
# 5 employees per farm verage from: https://www.nmpf.org/wp-content/uploads/immigration-survey-090915.pdf

PopSize = 100000 #Temp pop size for the purpose of testing the model

params = c(
  
  "beta.1"  = 0.000005,      #probability S15R0 -> infected by H1
  "beta.2"  = 0.00001,    #prob S15R0 -> infected by H5 (upon contact w sick cow)
  "beta.3"  = 0.001,    #prob I1R0 -> coinfected w H5 (upon contact w sick cow)
  "beta.4"  = 0.001,     #prob I5R0 -> coinfected w H1 
  "beta.5"  = 0.001,    #prob S5R1 -> infected by H5
  "beta.6"  = 0.001,      #prob S1R5 -> infected by H1
  # "beta.x"  = 0.0001,  #prob any infected by Hx (recombinant spillover)
  
  "vaxef"    = 0.1,    #vaccination effectiveness against H1 (%)
  "vaxrate"  = 0.05,    #rate of vaccine uptake (per day?)

  "gamma.1"  = 0.01,      #recovery rate from H1
  "gamma.5"  = 0.1,      #recovery rate from H5
  "gamma.2"  = 0.1,      #recovery rate from coinfection 
# "gamma.x"  = 1/10,     #recovery rate from Hx 
  
  "cow"      = 0         #no. sick cows
  
#  reassort  = 0
  
)


cowflu.model = function(t, x, params) {
  
  #State variables
  S15R0 = x[1]  #susceptible agricultural workers per county
  I1R0  = x[2]  #fully susceptible workers infected with H1
  I5R0  = x[3]  #fully susceptible workers infected with H5
  
  S5R1  = x[4]  #H1 resistant workers susceptible to H5
  I5R1  = x[5]  #H1 resistant workers infected with H5
  
  S1R5  = x[6]  #H5 resistant workers susceptible to H1
  I1R5  = x[7]  #H5 resistant workers infected with H1
  
  I15R0 = x[8]  #workers co-infected with both H1 and H5
  # IX    = x[9]  #workers infected with new spillover variant
  
  R15   = x[9] #fully resistant workers

  
  #Parameters (explained below)
  beta.1  = params["beta.1"]   #probability S15R0 -> infected by H1
  beta.2  = params["beta.2"]   #prob S15R0 -> infected by H5 (upon contact w sick cow)
  beta.3  = params["beta.3"]   #prob I1R0 -> coinfected w H5 (upon contact w sick cow)
  beta.4  = params["beta.4"]   #prob I5R0 -> coinfected w H1 
  beta.5  = params["beta.5"]   #prob S5R1 -> infected by H5 (upon contact w sick cow)
  beta.6  = params["beta.6"]   #prob S1R5 -> infected by H1
  # beta.x  = params["beta.x"]   #prob any infected by Hx (recombinant spillover)
   
  vaxef    = params["vaxef"]   #vaccination effectiveness against H1 (%)
  vaxrate  = params["vaxrate"]   #rate of vaccine uptake
  
  cow      = params["cow"]  #contact rate with sick cows
  
  gamma.1  = params["gamma.1"]  #recovery rate from H1
  gamma.5  = params["gamma.5"]  #recovery rate from H5
  gamma.2  = params["gamma.2"]  #recovery rate from coinfection 
  
  # reassort = params[14]  #probability H1 and H5 reassort into Hx (new variant)

  #System of equations
  
  dS15R0t <-   - beta.1*S15R0*I1R0 - S15R0*beta.1*I15R0 - S15R0*beta.1*I1R5 - S15R0*cow*beta.2 - S15R0*vaxrate*vaxef # H5 infection = bovine-to-human, no human-to-human transmission
  
  dI1R0t  <-   beta.1*S15R0*I1R0 + S15R0*beta.1*I15R0 + S15R0*beta.1*I1R5 - I1R0*gamma.1 - I1R0*cow*beta.3
  
  dI5R0t  <-   S15R0*cow*beta.2 - I5R0*beta.4*I1R0 - I5R0*beta.4*I15R0 - I5R0*beta.4*I1R5 - I5R0*gamma.5 #Problem child lol
  
  dS5R1t  <-   I1R0*gamma.1 - S5R1*cow*beta.5 + S15R0*vaxrate*vaxef 
  
  dI5R1t  <-   S5R1*cow*beta.5 - I5R1*gamma.5
  
  dS1R5t  <-   I5R0*gamma.5 - S1R5*beta.6*I1R0 - S1R5*beta.6*I15R0 - S1R5*beta.6*I1R5 
    
  dI1R5t  <-   S1R5*beta.6*I1R0 + S1R5*beta.6*I15R0 + S1R5*beta.6*I1R5 - I1R5*gamma.1 
  
  dI15R0t <-   I1R0*cow*beta.3 + I5R0*beta.4*I1R0 + I5R0*beta.4*I15R0 + I5R0*beta.4*I1R5 - I15R0*gamma.2
    
  dR15t   <-   I5R1*gamma.5 + I1R5*gamma.1 + I15R0*gamma.2
  
  
  output = c(dS15R0t, dI1R0t, dI5R0t, dS5R1t, dI5R1t, dS1R5t, dI1R5t, dI15R0t, dR15t)
  names(output) = c('S15R0', 'I1R0', 'I5R0', 'S5R1', 'I5R1', 'S1R5', 'I1R5', 'I15R0', 'R15')
  
  return(list(output))
}

# Initial conditions
x0        = numeric(9)

x0[1]     = PopSize
x0[2]     = 1
x0[3]     = 0
x0[3]     = 0
x0[4]     = 0
x0[5]     = 0
x0[6]     = 0
x0[7]     = 0
x0[8]     = 0
x0[9]     = 0
names(x0) = c('S15R0', 'I1R0', 'I5R0', 'S5R1', 'I5R1', 'S1R5', 'I1R5', 'I15R0', 'R15')

x0[]

# Run the model 
odeSim1 = ode(y = x0,0:100, cowflu.model, params)#, method = "radau")


# # Plot results
matplot(odeSim1[,1], odeSim1[,3:9])

# More detailed plot 
Simdat <- data.frame(odeSim1)

Simdat_long <-
  Simdat %>%
  pivot_longer(!time, names_to = "groups", values_to = "count")

ggplot(Simdat_long, aes(x= time, y = count, group = groups, colour = groups, )) + 
  geom_line(size = 1.2) 


# ####
# 
# reassort1 = 1/100000
# 
# prob = odeSim1[,7] * reassort1
# outbreak.indicator.vector = array(100000)
# 
# for (i in 1:100000){
# 
# draws = runif(length(prob))
# 
#   reassort.index = draws<prob
# 
#   outbreak.indicator = sum(reassort.index) > 0
# 
#   outbreak.indicator.vector[i] = outbreak.indicator
# 
#   }
# 
# 
# sum(outbreak.indicator.vector)/100000 #fraction of scenarios in which there is a recomb event
# 
# ####








#The simplest parts of the model (S15R0 and I1R0)

#What happens to the fully susceptible population?
# 
Suscept <- Simdat_long %>%
  filter(groups == "S15R0" | groups == "I1R0")

ggplot(Suscept, aes(x= time, y = count, group = groups, colour = groups )) +
  geom_line(size = 1)

# 
# #Focus on Coinfection
# HCo.plot <-
#   Simdat_long %>%
#   filter(groups == "I15R0")
# 
# ggplot(HCo.plot, aes(x= time, y = count, group = groups, colour = groups)) +
#   geom_line()
# 
# #Focus on Coinfection and Hx (new variant)
# Hx.plot <-
#   Simdat_long %>%
#   filter(groups == "IxR0" | groups == "I15R0")
# 
# ggplot(Hx.plot, aes(x= time, y = count, group = groups, colour = groups)) +
#   geom_line()
# 
# 
# #Focus on Infections (NO Coinfection)
# Inf.plot <-
# Simdat_long %>%
#   filter(groups == "I1R0" | groups == "I5R0" | groups == "I5R1" | groups == "I1R5" | groups == "IxR0")
# 
# ggplot(Inf.plot, aes(x= time, y = count, group = groups, colour = groups)) +
#   geom_line()
# 
# #Focus on H5 infections
# H5.inf <-
#   Simdat_long %>%
#   filter(groups == "I5R0" | groups == "I5R1")
# 
# ggplot(H5.inf, aes(x= time, y = count, group = groups, colour = groups)) +
#   geom_line()
# 
# 
# 
# 
