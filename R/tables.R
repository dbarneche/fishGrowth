tableRounding  <-  function (x, rounding = 2) {
    format(round(x, rounding), nsmall = rounding)
}

makeTable1  <-  function (geStanOut, gfStanOut) {
	geValues  <-  c('scalingAlpha_Intercept', 'Er', 'Ei_Intercept', 'Topt_Intercept', 'lnBoTs_Intercept', 'sd(lnBoTs_Intercept)')
	gfValues  <-  c('lnMopt', 'invKT', 'Intercept', 'sd(Intercept)', 'sd(lnMopt)', 'sd(invKT)', 'cor(Intercept,lnMopt)', 'cor(Intercept,invKT)', 'cor(lnMopt,invKT)')

	table1  <-  cbind(rbind(
							tableRounding(as.matrix(geStanOut[geValues[1:5], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2), 
							matrix(rep('', 6), 2, 3),
							tableRounding(as.matrix(geStanOut[geValues[6], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2), 
							matrix(rep('', 15), 5, 3)
							),
					  rbind(
							tableRounding(as.matrix(gfStanOut[gfValues[1:2], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2), 
							matrix(rep('', 6), 2, 3),
							tableRounding(as.matrix(gfStanOut[gfValues[3], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2), 
							matrix(rep('', 6), 2, 3),
							tableRounding(as.matrix(gfStanOut[gfValues[4:9], c('Estimate', 'l-95% CI', 'u-95% CI')]), 2)
							)
					  )
	rownames(table1)  <-  c(geValues[1:5], '', '', gfValues[4:9])
	table1
}

writeTable1  <-  function (dest, table1) {
	write.csv(table1, dest)
}

makeTableS1  <-  function (fullModel, simplerModel, simplestModel) {
	pFull      <-  c('b_scalingAlpha_Intercept', 'b_logitEr_Intercept', 'b_Ei_Intercept', 'b_Topt_Intercept', 'b_lnBoTs_Intercept', 'b_trophicLevelSlope_Intercept', 'b_aspectRatioSlope_Intercept', 'sd_family__lnBoTs_Intercept', 'sd_family__scalingAlpha_Intercept', 'sd_family__logitEr_Intercept', 'sd_family__Topt_Intercept', 'cor_family__lnBoTs_Intercept__scalingAlpha_Intercept', 'cor_family__lnBoTs_Intercept__logitEr_Intercept', 'cor_family__lnBoTs_Intercept__Topt_Intercept', 'cor_family__scalingAlpha_Intercept__logitEr_Intercept', 'cor_family__scalingAlpha_Intercept__Topt_Intercept', 'cor_family__logitEr_Intercept__Topt_Intercept')
	pSimpler   <-  pFull[-7]
	pSimplest  <-  pFull[-(6:7)]

    brmsOutFull      <-  t(brms::posterior_samples(fullModel, pars = pFull, exact_match = TRUE, as.matrix = TRUE))[pFull, ]
    brmsOutSimpler   <-  t(brms::posterior_samples(simplerModel, pars = pSimpler, exact_match = TRUE, as.matrix = TRUE))[pSimpler, ]
    brmsOutSimplest  <-  t(brms::posterior_samples(simplestModel, pars = pSimplest, exact_match = TRUE, as.matrix = TRUE))[pSimplest, ]
	
    # transform logitEr values to Er
    transformLogitEr  <-  function (mat) {
    	mat['b_logitEr_Intercept', ]  <-  mat['b_Ei_Intercept', ] / (1 + exp(-1 * mat['b_logitEr_Intercept', ]))
    	mat	
    }
	brmsOutFull      <-  transformLogitEr(brmsOutFull)
	brmsOutSimpler   <-  transformLogitEr(brmsOutSimpler)
	brmsOutSimplest  <-  transformLogitEr(brmsOutSimplest)

	sumTab  <-  function (output) {
		plyr::adply(output, 1, function (x) {
			qts  <-  quantile(x, probs = c(0.025, 0.975))
			data.frame('Estimate' = mean(x), 'l95_CI' = qts[1], 'u95_CI' = qts[2], stringsAsFactors =  FALSE)
		})
	}
	
	# transform intercepts to J / g
	transformBoToJ  <-  function (mat) {
		mat[mat$parameters == 'b_lnBoTs_Intercept', c('Estimate', 'l95_CI', 'u95_CI')]  <-  mat[mat$parameters == 'b_lnBoTs_Intercept', c('Estimate', 'l95_CI', 'u95_CI')] + log(39e3)
		mat
	}
	sumTabFull      <-  transformBoToJ(sumTab(brmsOutFull))
	sumTabSimpler   <-  transformBoToJ(sumTab(brmsOutSimpler))
	sumTabSimplest  <-  transformBoToJ(sumTab(brmsOutSimplest))

	cbind(
		  rbind(tableRounding(as.matrix(sumTabFull[1:7, 2:4]), 2), matrix(rep('', 6), 2, 3), tableRounding(as.matrix(sumTabFull[8:nrow(sumTabFull), 2:4]), 2)),
		  rbind(tableRounding(as.matrix(sumTabSimpler[1:6, 2:4]), 2), matrix(c(rep('-', 3), rep('', 6)), 3, 3, byrow = TRUE), tableRounding(as.matrix(sumTabSimpler[7:nrow(sumTabSimpler), 2:4]), 2)),
		  rbind(tableRounding(as.matrix(sumTabSimplest[1:5, 2:4]), 2), matrix(c(rep('-', 6), rep('', 6)), 4, 3, byrow = TRUE), tableRounding(as.matrix(sumTabSimplest[6:nrow(sumTabSimplest), 2:4]), 2))
		  )
}

writeTableS1  <-  function (dest, tableS1) {
	write.csv(tableS1, dest)
}

writeTableS2  <-  function (dest) {
	tableS2  <-  data.frame('Original units' = c('joules / d', 'nL O2 / h', 'uL O2 / h', 'mL O2 / h', 'ug O2 / h', 'nmol O2 / h', 'umol O2 / h', 'mg O2 / d'), 'Multiplication factor to yield mg O2 / h' = c(0.00284616, 0.000001429, 0.001429, 1.429, 0.001, 0.000032, 0.032, 0.041666667), stringsAsFactors = FALSE, check.names = FALSE)
	write.csv(tableS2, dest, row.names = FALSE)
}
