rm(list=ls())
source('R/functions-analyses.R')

oozeki       <-  readFile('data/experimentalData/oozeki_y_1994_j41_207_fig1.csv')
oozeki$days  <-  round(oozeki$days)
avMetRate    <-  tapply(oozeki$ulO2_h, oozeki$days, mean)
data         <-  data.frame(days=as.numeric(names(avMetRate)), ulO2_h=avMetRate, stringsAsFactors=FALSE)

# according to oozeki_y_1992_j39_59, the equation describing dry mass (mg) is
# dry mass (mg) = 0.00391 * 10^(0.105 * days)

data$dry_mass_mg  <-  0.00391 * 10^(0.105 * data$days)
data$growth_mg_d  <-  c(NA, instantaneousGrowth(data$dry_mass_mg[2:nrow(data)],
							   		            data$dry_mass_mg[1:(nrow(data)-1)],
									            data$days[2:nrow(data)],
									            data$days[1:(nrow(data)-1)]))
