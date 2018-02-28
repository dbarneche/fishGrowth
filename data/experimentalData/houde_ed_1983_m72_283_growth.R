#########################################
# ALL CONTAINED IN TABLES 3-5, PAGE 288 
# MASS IS GIVEN IN DRY WEIGHT (µg)
# GROWTH IS GIVEN IN CAL/DAY
# CONVERSION: 0.0050 CAL/µg OF DRY WEIGHT
# TEMPERATURES ARE ASSUMED TO BE THE SAME 
# AS REPORTED IN TABLE 2, PAGE 287
# BECAUSE THE AUTHORS USE A MODEL FROM
# ANOTHER OF THEIR OWN PAPERS, I 
# JUST REAPPLIED THE MODEL TO ESTIMATE
# THE PARAMETERS AND RECALCULATE THE 
# GROWTH FOR THE MASS VALUES FOR WHICH WE
# HAVE METABOLIC RATE MEASUREMENTS
#########################################

# ANCHOA MITCHILLI, TABLE 3
G  <-  c(0.021, 0.031, 0.041, 0.103, 0.206) / 0.005 # transform growth to µg/day
W  <-  c(10,15,20,50,100) # mass reported on table
M  <-  coef(summary(lm(G ~ W)))
N  <-  c(8.9, 31.3, 98.5, 424.4) # actual mass for which metabolic rates were measured 

M['(Intercept)','Estimate'] + M['W','Estimate'] * N 
3.689726
12.910395
40.572402
174.724904

# ARCHOSARGUS RHOMBOIDALIS, TABLE 4
G  <-  c(0.025, 0.037, 0.050, 0.125, 0.250) / 0.005 # transform growth to µg/day
W  <-  c(10,15,20,50,100) # mass reported on table
M  <-  coef(summary(lm(G ~ W)))
N  <-  c(18.1, 41.9, 66.2) # actual mass for which metabolic rates were measured 

M['(Intercept)','Estimate'] + M['W','Estimate'] * N 
9.021075 
20.931238 
33.091616

# ACHIRUS LINEATUS, TABLE 5
G  <-  c(0.019, 0.029, 0.038, 0.096, 0.193) / 0.005 # transform growth to µg/day
W  <-  c(10,15,20,50,100) # mass reported on table
M  <-  coef(summary(lm(G ~ W)))
N  <-  c(14.3, 18.1, 23.1, 49.4, 63.5, 248.4) # actual mass for which metabolic rates were measured 

M['(Intercept)','Estimate'] + M['W','Estimate'] * N 
5.454021
6.922633
8.855018
19.019359
24.468683
95.928256
