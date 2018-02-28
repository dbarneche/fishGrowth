library(plyr)

rm(list=ls())
source('R/functions-analyses.R')

growth  <-  readFile('data/experimentalData/rombough_pj_1994_f8_178_fig1.csv')
rates   <-  readFile('data/experimentalData/rombough_pj_1994_f8_178_fig3.csv')

growth$days    <-  ceiling(growth$days)
rates$days     <-  ceiling(rates$days)
rates$mass_mg  <-  NA

temps  <-  unique(growth$temp_celsius)
n      <-  0

for(i in temps) {
	g1  <-  growth[growth$temp_celsius == i, ]	
	r1  <-  rates[rates$temp_celsius == i, ]

	for(k in 1:nrow(r1)) {
		n       <-  n + 1
		diffe   <-  abs(g1$days - r1$days[k])
		j       <-  which.min(diffe)
		if(diffe[j] <= 1) {
			rates$mass_mg[n]  <-  g1$mass_mg[j]
		}

	}
}

rates  <-  rates[complete.cases(rates), ]

instantaneousGrowth  <-  function(m1, m0, t1, t0) {
	(log(m1/m0) * m1)/ (t1 - t0)
}

calcGrowth  <-  function(data) {
	data$growth_mg_day  <-  NA
	for(j in 2:nrow(data)) {
		data$growth_mg_day[j]  <-  instantaneousGrowth(data$mass_mg[j], data$mass_mg[j-1], data$days[j], data$days[j-1])
	}
	data
}

rates  <-  ddply(rates, .(temp_celsius), calcGrowth)
rates  <-  rates[complete.cases(rates), ]

write.csv(rates, 'data/experimentalData/rombough_pj_1994_f8_178_fig1_3_merged.csv', row.names=FALSE)

rates           <-  rates[rates$growth_mg_day > 0, ]
rates$lnRate    <-  log(rates$ug_o2_h)
rates$lnGrowth  <-  log(rates$growth_mg_day)
rates$lnMass    <-  log(rates$mass_mg)
rates$invKT     <-  1/8.62e-5*(1/281.15 - 1/(rates$temp_celsius + 273.15))

lmG    <-  lm(lnGrowth ~ lnMass + invKT, data=rates)
lmR    <-  lm(lnRate ~ lnMass + invKT, data=rates)

massCorrectedRate  <-  rates$lnRate - coef(lmR)[2]*rates$lnMass
tempCorrectedRate  <-  rates$lnRate - coef(lmR)[3]*rates$invKT
plot(massCorrectedRate ~ rates$invKT)
plot(tempCorrectedRate ~ rates$lnMass)

massCorrectedGrowth  <-  rates$lnGrowth - coef(lmG)[2]*rates$lnMass
tempCorrectedGrowth  <-  rates$lnGrowth - coef(lmG)[3]*rates$invKT
plot(massCorrectedGrowth ~ rates$invKT)
plot(tempCorrectedGrowth ~ rates$lnMass)

plot(lnGrowth ~ lnRate, data=rates)

