######################
# AUXILLIARY FUNCTIONS
######################
extrafont::loadfonts(quiet = TRUE)

toDev  <-  function (expr, dev, filename, ..., verbose=TRUE) {
    if ( verbose ) {
        cat(sprintf('Creating %s\n', filename))
    }
    dev(filename, family = 'CM Roman', ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}

toPdf  <-  function (expr, filename, ...) {
    toDev(expr, pdf, filename, ...)
}

mainLabel  <-  function (intercept, slope, ln = FALSE, negative = FALSE, ...) {
    if (negative) {
        if (ln) {
            substitute('mean trend: ln('~italic(y)~')' == B - A * italic(x), list(B = LoLinR::rounded(intercept, ...), A=LoLinR::rounded(-1 * slope, ...)))
        } else {
            substitute('mean trend: '~italic(y) == B - A * italic(x), list(B = LoLinR::rounded(intercept, ...), A = LoLinR::rounded(-1 * slope, ...)))
        }
    } else {
        if (ln) {
            substitute('mean trend: ln('~italic(y)~')' == B + A * italic(x), list(B = LoLinR::rounded(intercept, ...), A = LoLinR::rounded(slope, ...)))
        } else {
            substitute('mean trend: '~italic(y) == B + A * italic(x), list(B = LoLinR::rounded(intercept, ...), A = LoLinR::rounded(slope, ...)))
        }
    }
}

ogm  <-  function (mass_time, Bo, energiesVector, alpha, asympMass) {
    # mass as a function of mass_time over ontogeny starting at m_o = 0
    # i.e. m_o / asympMass = 0 in eqn S3
    exponent  <-  -1 * (Bo / energiesVector$Em) * (1 - alpha) * mass_time * asympMass^(alpha - 1)
    asympMass * (1 - exp(exponent))^(1 / (1 - alpha))
}

ogmTime  <-  function (mass, Bo, energiesVector, alpha, asympMass) {
    # time required to reach mass starting from mass mass = 0
    # this is the function ogm above (eqn S3) solved for time t
    -1 * (energiesVector$Em * log(1 - (mass / asympMass)^(1 - alpha))) / (Bo * (1 - alpha) * asympMass^(alpha - 1))
}

assimilate  <-  function (mass, Bo, fRatio, energiesVector, alpha, asympMass) {
    # rate of assimilation at mass m, eqn 5
    Bo * mass^alpha * (fRatio + (energiesVector$Ec / energiesVector$Em) * (1 - (mass / asympMass)^(1 - alpha)))
}

efficiencyFunction  <-  function (mass_time, Bo, fRatio, energiesVector, alpha, asympMass) {
    # wrapper function
    ontMass  <-  ogm(mass_time, Bo, energiesVector, alpha, asympMass)
    assimilate(ontMass, Bo, fRatio, energiesVector, alpha, asympMass)
}

efficiency  <-  function (mass, Bo, fRatio, energiesVector, alpha, asympMass) {
    # eqn 6
    epsilon  <-  rep(0, length(mass))
    for (i in seq_along(epsilon)) {
        ontTime     <-  ogmTime(mass = mass[i], Bo = Bo, energiesVector = energiesVector, alpha = alpha, asympMass = asympMass)
        epsilon[i]  <-  energiesVector$Ec * mass[i] / integrate(efficiencyFunction, lower = 0, upper = ontTime, Bo = Bo, fRatio = fRatio, energiesVector = energiesVector, alpha = alpha, asympMass = asympMass)$value
    }
    epsilon
}

getTe  <-  function (exponent, ppmr) {
    exp((exponent - 0.25) * log(ppmr))
}

trianglePolygon  <-  function (xMin, xMax, yMin, yMax, cs, invert = FALSE) {
    # first figure out slope
    xBar      <-  xMin + (xMax - xMin) / 2
    if (invert) {
        slope     <-  (yMin - yMax) / (xMax - xBar)
        intcp     <-  yMax - slope * xBar
        brks      <-  seq(yMax, yMin, length.out = 6)
        xs        <-  (brks - intcp) / slope
        interval  <-  (xs[2] - xBar)
    } else {
        slope     <-  (yMax - yMin) / (xBar - xMin)
        intcp     <-  yMin - slope * xMin
        brks      <-  seq(yMin, yMax, length.out = 6)
        xs        <-  (brks - intcp) / slope
        interval  <-  (xs[2] - xMin)
    }
    
    #base
    polygon(c(xMin, xMax, xMax - interval, xMin + interval, xMin), c(brks[1], brks[1], brks[2], brks[2], brks[1]), col = LoLinR::transparentColor(cs, 0.5), border = NA)
    #middle
    xMinM  <-  xMin + interval * 2
    xMaxM  <-  xMax - interval * 2
    polygon(c(xMinM, xMaxM, xMaxM - interval, xMinM+interval, xMinM), c(brks[3], brks[3], brks[4], brks[4], brks[3]), col = LoLinR::transparentColor(cs, 0.5), border = NA)
    #top
    xMinT  <-  xMin + interval * 4
    xMaxT  <-  xMax - interval * 4
    polygon(c(xMinT, xMaxT, xBar, xMinT), c(brks[5], brks[5], brks[6], brks[5]), col = LoLinR::transparentColor(cs, 0.5), border = NA)
}

RGBcomponents  <- function (myColor) {
    col2rgb(myColor)[, , drop = TRUE]
}

rgbDistance  <-  function (rgbDiffs) {
    sqrt(rgbDiffs[1]^2 + rgbDiffs[2]^2 + rgbDiffs[3]^2)
}

plotPreset  <-  function () {
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
}

###############
# PAPER FIGURES
###############
makeFigure1  <-  function (dest, ...) {
    toPdf(fig1(...), dest, width = 11, height = 6.5)
    extrafont::embed_fonts(dest)
}

fig1  <-  function (growthFishBase, gfStanOut, growthEmData, geStanOut) {
    ##########################
    # FISHBASE DATA CORRECTION
    ##########################
    alphaFixedEffectsFishBase       <-  gfStanOut['lnMopt', 'Estimate']
    erFixedEffectsFishBase          <-  gfStanOut['invKT', 'Estimate']
    alphaRandomEffectsFishBase      <-  gfStanOut[paste0('lnMopt_', growthFishBase$family), 'Estimate']
    erRandomEffectsFishBase         <-  gfStanOut[paste0('invKT_', growthFishBase$family), 'Estimate']
    alphaVectorFishBase             <-  alphaFixedEffectsFishBase + alphaRandomEffectsFishBase
    erVectorFishBase                <-  erFixedEffectsFishBase + erRandomEffectsFishBase
    growthFishBase$correctedlnRate  <-  growthFishBase$lnMaxGrowths - erVectorFishBase * growthFishBase$invKT
    goptCorrFactor                  <-  goptCorrection(0.77)
    growthFishBase$correctedlnRate  <-  log(exp(growthFishBase$correctedlnRate) / goptCorrFactor)
    
    #####################
    # LAB DATA CORRECTION
    #####################
    experimentalGrowth                  <-  growthEmData
    experimentalGrowth$sppNum           <-  as.numeric(as.factor(experimentalGrowth$species))
    alphaFixedEffectsExperimental       <-  geStanOut['scalingAlpha_Intercept', 'Estimate']
    logitErFixedEffectsExperimental     <-  geStanOut['logitEr', 'Estimate']
    erFixedEffectsExperimental          <-  geStanOut['Er', 'Estimate']
    eiFixedEffectsExperimental          <-  geStanOut['Ei', 'Estimate']
    toptFixedEffectsExperimental        <-  geStanOut['Topt', 'Estimate']
    experimentalGrowth$correctedlnRate  <-  experimentalGrowth$lnRate - erFixedEffectsExperimental / 8.62e-5 * (1/288.15 - 1/experimentalGrowth$growthTempKelvin) - log(1 + erFixedEffectsExperimental/(eiFixedEffectsExperimental - erFixedEffectsExperimental) * exp(eiFixedEffectsExperimental / 8.62e-5 * (1 / toptFixedEffectsExperimental - 1 / experimentalGrowth$growthTempKelvin)))
    
    lnMass          <-  c(growthFishBase$lnMopt, experimentalGrowth$lnMass)
    lnGrowth        <-  c(growthFishBase$correctedlnRate, experimentalGrowth$correctedlnRate)
    cols            <-  c(rep('grey90', nrow(growthFishBase)), rep('white', nrow(experimentalGrowth)))
    pchs            <-  c(rep(21, nrow(growthFishBase)), rep(22, nrow(experimentalGrowth)))
    
    par(mfrow = c(1, 2), mai = c(1.2, 0.82, 0.82, 0.60), omi = c(1, 1, 1, 1), cex = 1, cex.lab = 1.3, mgp = c(3, 0.5, 0), tck = -0.03)
    plot(NA, xlab = '', ylab = '', xlim = c(-10, 14), ylim = c(-13, 8), axes = FALSE)
    axis(1, at = seq(-10, 14, 8), cex.axis = 1.1)
    axis(2, las = 1, cex.axis = 1.1)
    axis(3, at = seq(-10, 14, 8), labels = c(expression(paste('4.5', ' x ', 10^-5, sep = '')),
                                             expression(paste('1.4', ' x ', 10^-2, sep = '')),
                                             expression(paste('4.0', ' x ', 10^2, sep = '')),
                                             expression(paste('1.2', ' x ', 10^6, sep = ''))))
    box()
    points(lnGrowth ~ lnMass, pch = pchs, col = 'grey60', bg = cols, cex = 0.9, lwd = 0.5)
    LoLinR::proportionalLabel(-0.2, 0.5, expression(paste('ln Growth rate @ 15 ' * degree * 'C (g d'^{-1}, ')'), sep = ''), adj = c(0.5, 0.5), xpd = NA, srt = 90, cex = 1.3)
    LoLinR::proportionalLabel(0.5, -0.25, 'ln Mass (g)', adj = c(0.5, 0.5), xpd = NA, cex = 1.3)
    LoLinR::proportionalLabel(px = 0.5, py = 1.25, 'Mass (g)', xpd = NA, cex = 1.3, adj = c(0.5, 0.5))

    for (j in unique(growthFishBase$family)) {
        xpoints  <-  range(growthFishBase$lnMopt[growthFishBase$family == j])
        nresA    <-  gfStanOut['lnMopt', 'Estimate'] + gfStanOut[paste0('lnMopt_', j), 'Estimate']
        nresboTc <-  log(exp(gfStanOut['Intercept', 'Estimate'] + gfStanOut[paste0('Intercept_', j), 'Estimate']) / goptCorrFactor)
        expr     <-  nresboTc + nresA * xpoints
        lines(xpoints, expr, lty = 2, col = 'black', lwd = 1)
    }
    lines(range(growthFishBase$lnMopt), log(exp(gfStanOut['Intercept', 'Estimate']) / goptCorrFactor) + alphaFixedEffectsFishBase * range(growthFishBase$lnMopt), lwd = 1.5, col = 'black')
    
    for (j in unique(experimentalGrowth$species)) {
        xpoints  <-  range(experimentalGrowth$lnMass[experimentalGrowth$species == j])
        nresboTc <-  geStanOut['lnBoTs_Intercept', 'Estimate'] + geStanOut[paste0('lnBoTs_Intercept_', j), 'Estimate']
        expr     <-  nresboTc + alphaFixedEffectsExperimental * xpoints
        lines(xpoints, expr, lty = 2, col = 'black', lwd = 1)
    }
    lines(range(experimentalGrowth$lnMass), geStanOut['lnBoTs_Intercept', 'Estimate'] + alphaFixedEffectsExperimental * range(experimentalGrowth$lnMass), lwd = 1.5, col = 'black')
    
    LoLinR::proportionalLabel(px = c(0.82, 0.95), py = c(0.2, 0.2), text = FALSE, type = 'l', lty = 1, lwd = 1.5, col = 'black')
    LoLinR::proportionalLabel(px = c(0.82, 0.95), py = c(0.1, 0.1), text = FALSE, type = 'l', lty = 2, lwd = 1, col = 'black')
    LoLinR::proportionalLabel(px = 0.8, py = 0.2, 'mean trend', adj = c(1, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(px = 0.8, py = 0.1, 'taxon-level variation', adj = c(1, 0.5), cex = 0.9)
    
    # point legends 
    LoLinR::proportionalLabel(px = 0.03, py = 0.95, text = FALSE, col = 'grey60', bg = 'white', adj = c(1, 0.5), pch = 22, lwd = 0.5, cex = 0.9, xpd = NA)
    LoLinR::proportionalLabel(px = 0.03, py = 0.85, text = FALSE, col = 'grey60', bg = 'grey90', adj = c(1, 0.5), pch = 21, lwd = 0.5, cex = 0.9, xpd = NA)
    # text legends
    LoLinR::proportionalLabel(px = 0.06, py = 0.95, lab = substitute('I: Embryo/Larval rate (' * italic('n') * ' = ' * a * ')', list(a = nrow(experimentalGrowth))) , cex = 0.9, adj = c(0, 0.5), xpd = NA)
    LoLinR::proportionalLabel(px = 0.06, py = 0.85, lab = substitute('II: Max. rate, ' * italic('G'['opt']) * ' (' * italic('n') * ' = ' * a * ')', list(a = niceThousands(nrow(growthFishBase)))), cex = 0.9, adj = c(0, 0.5), xpd = NA)

    #####################
    # TEMPERATURE EFFECTS
    #####################
    growthFishBase$tempCorrectedlnRate      <-  growthFishBase$lnMaxGrowths - alphaVectorFishBase * growthFishBase$lnMopt
    growthFishBase$tempCorrectedlnRate      <-  log(exp(growthFishBase$tempCorrectedlnRate) / goptCorrFactor)
    experimentalGrowth$tempCorrectedlnRate  <-  experimentalGrowth$lnRate - alphaFixedEffectsExperimental * experimentalGrowth$lnMass
    
    invKT           <-  c(growthFishBase$invKT, experimentalGrowth$invKT)
    templnGrowth    <-  c(growthFishBase$tempCorrectedlnRate, experimentalGrowth$tempCorrectedlnRate)
    
    plot(NA, xlab = '', ylab = '', xlim = (1 / 8.62e-5 * (1 / 288.15 - 1 / (273.15 + c(-1, 30)))), ylim = c(-6, -1), axes = FALSE)
    axis(1, at = round(1 / 8.62e-5 * (1 / 288.15 - 1 / (273.15 + seq(0, 30, 10))), 2), cex.axis = 1.1)
    axis(2, las = 1, cex.axis = 1.1) 
    axis(3, at = round(1 / 8.62e-5 * (1 / 288.15 - 1 / (273.15 + seq(0, 30, 10))), 2), labels = seq(0, 30, 10))
    box()
    points(templnGrowth ~ invKT, pch = pchs, col = 'grey60', bg = cols, cex = 0.9, lwd = 0.5)
    LoLinR::proportionalLabel(-0.2, 0.5, expression(paste('ln Growth rate @ 1 g (g d'^{-1}, ')'), sep = ''), adj = c(0.5, 0.5), xpd = NA, srt = 90, cex = 1.3)
    LoLinR::proportionalLabel(0.5, -0.25, expression(paste('Inverse Temperature, 1/', italic(kT[s]), ' - 1/', italic(kT), ' (eV'^{-1}, ')', sep = '')), adj = c(0.5, 0.5), xpd = NA, cex = 1.3)
    LoLinR::proportionalLabel(px = 0.5, py = 1.25, expression(paste('Temperature (' * degree, 'C)', sep = '')), xpd = NA, cex = 1.3, adj = c(0.5, 0.5))
    
    for (j in unique(growthFishBase$family)) {
        xpoints  <-  range(growthFishBase$invKT[growthFishBase$family == j])
        nresEr   <-  gfStanOut['invKT', 'Estimate'] + gfStanOut[paste0('invKT_', j), 'Estimate']
        nresboTc <-  log(exp(gfStanOut['Intercept', 'Estimate'] + gfStanOut[paste0('Intercept_', j), 'Estimate']) / goptCorrFactor)
        expr     <-  nresboTc + nresEr * xpoints
        lines(xpoints, expr, lty = 2, col = 'black', lwd = 1)
    }
    lines(range(growthFishBase$invKT), log(exp(gfStanOut['Intercept', 'Estimate']) / goptCorrFactor) + erFixedEffectsFishBase * range(growthFishBase$invKT), lwd = 1.5, col = 'black')
    
    for (j in unique(experimentalGrowth$species)) {
        pts      <-  range(experimentalGrowth$growthTempKelvin)
        xpoints  <-  seq(pts[1], pts[2], length.out = 50)
        expr     <-  (geStanOut['lnBoTs_Intercept', 'Estimate']  + geStanOut[paste0('lnBoTs_Intercept_', j), 'Estimate']) + 
            erFixedEffectsExperimental / 8.62e-5 * (1 / 288.15 - 1 / xpoints) - log(1 + exp(eiFixedEffectsExperimental / 8.62e-5 * (1 / toptFixedEffectsExperimental - 1 / xpoints)) * erFixedEffectsExperimental / (eiFixedEffectsExperimental - erFixedEffectsExperimental))
        lines(1 / 8.62e-5 * (1 / 288.15 - 1 / xpoints), expr, lty = 2, col = 'black', lwd = 1)
    }
    meanNx         <-  1 / 8.62e-5 * (1 / 288.15 - 1 / min(experimentalGrowth$growthTempKelvin):max(experimentalGrowth$growthTempKelvin))
    meanExprNlmer  <-  geStanOut['lnBoTs_Intercept', 'Estimate'] + erFixedEffectsExperimental / 8.62e-5 * (1 / 288.15 - 1 / min(experimentalGrowth$growthTempKelvin):max(experimentalGrowth$growthTempKelvin)) - log(1 + exp(eiFixedEffectsExperimental / 8.62e-5 * (1 / toptFixedEffectsExperimental - 1 / min(experimentalGrowth$growthTempKelvin):max(experimentalGrowth$growthTempKelvin))) * erFixedEffectsExperimental / (eiFixedEffectsExperimental - erFixedEffectsExperimental))
    lines(meanNx, meanExprNlmer, col = 'black', lwd = 1.5)
}

makeFigure2  <-  function (dest, ...) {
    toPdf(fig2(...), dest, width = 7, height = 6.5)
    extrafont::embed_fonts(dest)
}

fig2  <-  function (fishBaseGrowthDataWithEm, growthEmData, tempDependenceEmFishBaseData, tempDependenceEmExperimentalData) {
    modelFishBase  <-  tempDependenceEmFishBaseData
    modelExpData   <-  tempDependenceEmExperimentalData

    fixefsFishBase  <-  brms::fixef(modelFishBase, pars = 'b')
    ranefsFishBase  <-  brms::ranef(modelFishBase, pars = 'b')
    fixefsExperim   <-  brms::fixef(modelExpData, pars = 'b')
    ranefsExperim   <-  brms::ranef(modelExpData, pars = 'b')

    fishBaseGrowthDataWithEm$correctedLnEm  <-  log(fishBaseGrowthDataWithEm$Em_J_g) - ranefsFishBase$family[fishBaseGrowthDataWithEm$family, 'Estimate', 'Intercept'] - ranefsFishBase$family[fishBaseGrowthDataWithEm$family, 'Estimate', 'temperature'] * fishBaseGrowthDataWithEm$temperature
    growthEmData$correctedLnEm  <-  log(growthEmData$emStar) - ranefsExperim$species[growthEmData$species, 'Estimate', 'Intercept'] - ranefsExperim$species[growthEmData$species, 'Estimate', 'tempC'] * growthEmData$tempC

    par(cex = 1, omi = c(1, 1, 0.5, 1), cex.axis = 1.2, cex.lab = 1.4, xpd = NA)
    plot(NA, xlim = c(-3, 33), ylim = c(4, 12.2), xlab = expression(paste('Temperature ('*degree,'C)', sep = '')), ylab = expression(paste('ln ', italic(E['m'])*' (J g'^{-1}, ')', sep = '')), axes = FALSE)
    axis(1)
    axis(2, las=1)
    box()

    points(correctedLnEm ~ temperature, data = fishBaseGrowthDataWithEm, pch = 21, col = 'grey60', bg = 'grey90', cex = 1.1, lwd = 0.5)
    points(correctedLnEm ~ tempC, data = growthEmData, pch = 22, col = 'grey60', bg = 'white', cex = 1.1, lwd = 0.5)

    LoLinR::proportionalLabel(0.03, 0.95, text = FALSE, pch = 22, col = 'grey60', bg = 'white', cex = 1.1, lwd = 0.5)
    LoLinR::proportionalLabel(0.03, 0.85, text = FALSE, pch = 21, col = 'grey60', bg = 'grey90', cex = 1.1, lwd = 0.5)
    LoLinR::proportionalLabel(0.05, 0.95, substitute(z * italic(y) == b + italic(x) %.% a * ' (N.S.)', list(z = substitute(italic('n') == k * ' obs; ', list(k = nrow(growthEmData))), a = LoLinR::rounded(fixefsExperim['tempC', 'Estimate'], 2), b = LoLinR::rounded(fixefsExperim['Intercept', 'Estimate'], 2))), adj = c(0, 0.5))
    LoLinR::proportionalLabel(0.05, 0.85, substitute(z * italic(y) == b + italic(x) %.% a, list(z = substitute(italic('n') == k * ' obs; ', list(k = nrow(fishBaseGrowthDataWithEm))), a = LoLinR::rounded(fixefsFishBase['temperature', 'Estimate'], 2), b = LoLinR::rounded(fixefsFishBase['Intercept', 'Estimate'], 2))), adj = c(0, 0.5))

}

makeFigure3  <-  function (dest) {
    toPdf(fig3(), dest, width = 7, height = 6.5)
    extrafont::embed_fonts(dest)
}

fig3  <-  function () {
    Bo              <-  exp(-6.23) * 2880e3 / 6 / 12 # Bo for standard rates
    alpha           <-  0.75
    asympMass       <-  100 # arbitrary asymptotic size
    mass            <-  seq(1, 98, length.out = 600) # arbitrary mass, so it's easy to calculate mass / asympMass
    energiesVector  <-  list(Em = NA, Ec = 24000 * 0.15) # Ec value takes into account water content of biomass 
    
    par(cex = 1, omi = c(1, 1, 0.5, 1), cex.axis = 1.2, cex.lab = 1.4, xpd = NA)
    plot(NA, ylim = c(0.005, 2), xlim = c(0, 1), xlab = expression(paste('Ontogenetic stage (', italic('m'), '/', italic('M'), ')', sep = '')), ylab = expression(paste('Transfer efficiency, ' * epsilon, sep = '')), axes = FALSE, log = 'y')
    usr  <-  par('usr')
    rect(usr[1], 10^usr[3], usr[2], 10^usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid(log = 'y')
    box()
    axis(1)
    axis(2, at = c(0.005, 0.02, 0.05, 0.2, 0.5, 1), labels = c('0.005', '0.02', '0.05', '0.2', '0.5', '1'), las = 1)
    
    x    <-  mass / 100
    Ems  <-  c(1e3, 5e3, 1e4)
    # f lies between 1 (i.e. no activity, b_active = 0) and 4
    # here we explore these boundaries
    for (i in seq_along(Ems)) {
        energiesVector$Em  <-  Ems[i]
        upperLine  <-  efficiency(mass = mass, Bo = Bo, fRatio = 1, energiesVector = energiesVector, alpha = alpha, asympMass = asympMass)
        lowerLine  <-  efficiency(mass = mass, Bo = Bo, fRatio = 4, energiesVector = energiesVector, alpha = alpha, asympMass = asympMass)
        polygon(c(x, rev(x), x[1]), c(lowerLine, rev(upperLine), lowerLine[1]), lty = i)
    }

    LoLinR::proportionalLabel(0.95, 0.95, substitute(italic('E'['m']) * ' (kJ g'^-1 * '):'), adj = c(1, 0.5), log = 'y')
    LoLinR::proportionalLabel(0.95, 0.88, '1', adj = c(0.5, 0.5), log = 'y')
    LoLinR::proportionalLabel(0.95, 0.81, '5', adj = c(0.5, 0.5), log = 'y')
    LoLinR::proportionalLabel(0.95, 0.74, '10', adj = c(0.5, 0.5), log = 'y')
    LoLinR::proportionalLabel(c(0.82, 0.91), c(0.88, 0.88), text = FALSE, type = 'l', lty = 1, log = 'y')
    LoLinR::proportionalLabel(c(0.82, 0.91), c(0.81, 0.81), text = FALSE, type = 'l', lty = 2, log = 'y')
    LoLinR::proportionalLabel(c(0.82, 0.91), c(0.74, 0.74), text = FALSE, type = 'l', lty = 3, log = 'y')
}

makeFigure4  <-  function (dest) {
    toPdf(fig4(), dest, width = 7, height = 7.5)
    extrafont::embed_fonts(dest)
}

fig4  <-  function () {
    # basic parameter values
    Bo              <-  exp(-6.23) * 2880e3 / 6 / 12 # Bo for standard rates
    energiesVector  <-  list(Em = NA, Ec = 24000 * 0.15) # Ec value takes into account water content of biomass
    alpha           <-  0.75
    asympMass       <-  100  # arbitrary asymptotic size    
    fRatio          <-  2.4

    # varying values for looping
    Ems      <-  seq(1, 7e3, length.out = 600)
    mass     <-  seq(1, 71, length.out = 600) # arbitrary mass, so it's easy to calculate mass / asympMass
    allCols  <-  vector()
    
    par(cex = 1, omi = c(1, 1, 1.5, 1), cex.axis = 1.2, cex.lab = 1.4, xpd = NA)
    plot(NA, xlim = c(1, 7e3), ylim = c(1, 71), ylab = expression(paste('Ontogenetic stage (', italic('m'), '/', italic('M'), ')', sep = '')), xlab = expression(paste(italic(E['m']),' (J g'^{-1}, ')', sep = '')), xaxs = 'i', yaxs = 'i', axes = FALSE, cex.lab = 1.4, xpd = NA)
    
    effs     <-  data.frame()
    for (j in seq_along(Ems)) {
        energiesVector$Em  <-  Ems[j]
        effs   <-  rbind(effs, data.frame(Em = Ems[j], mass = mass, efficiency = round(efficiency(mass = mass, Bo = Bo, fRatio = fRatio, energiesVector = energiesVector, alpha = alpha, asympMass = asympMass), 3)))
    }
    allEffs  <-  round(seq(min(effs$efficiency), max(effs$efficiency), 0.001), 3)
    cols     <-  rev(colorRampPalette(c(colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Reds')[3:8]))(22), colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues')[5:7])(3)))(length(allEffs)))
    newCols  <-  cols[match(effs$efficiency, allEffs)]
    points(effs$Em, effs$mass, col = newCols, cex = 0.1, pch = 15)
    expos    <-  c(-0.1, 0, 0.1)
    tippEff  <-  vector(mode = 'numeric', length = length(expos))
    ltys     <-  c(2, 1, 3)
    for (k in seq_along(expos)) {
        tippEff[k]  <-  round(getTe(expos[k], 2327), 3)
        desired     <-  effs[effs$efficiency == tippEff[k], ]
        draw        <-  plyr::ddply(desired, .(mass), function (x) x[which.min(x$Em)[1], ])
        lines(draw$Em, draw$mass, lty = ltys[k])
    }
    
    xlbs  <-  c(0.01, 0.25, 0.62, 0.82)
    ylbs  <-  c(0.05, 0.4, 0.6, 0.9)
    exps  <-  c('>0.1', '>0', '<0', '<-0.1')
    for (i in seq_along(xlbs)) {
        LoLinR::proportionalLabel(xlbs[i], ylbs[i], substitute(italic('W') %prop% italic('m')^{z}, list(z = exps[i])), adj = c(0, 0.5), cex = 0.75)
    }
    trianglePolygon(xMin = 300, xMax = 2300, yMin = 5, yMax = 20, 'tomato4', invert = TRUE)
    trianglePolygon(xMin = 3300, xMax = 5300, yMin = 48, yMax = 63, 'dodgerblue4')
    
    axis(1, at = c(500, 2500, 4500, 6500))
    axis(2, at = seq(10, 70, 20), label = seq(0.1, 0.7, 0.2), las = 1)
    box()
    
    eqCols  <-  rev(colorRampPalette(c(colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Reds')[3:8]))(22), colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues')[5:7])(3)))(length(allEffs)))
    
    eqCols  <-  eqCols[allEffs <= tippEff[2] + (tippEff[2] - min(allEffs))]
    lseq    <-  seq(0.25, 0.75, length.out = length(eqCols))
    for (i in seq_along(lseq)) {
        LoLinR::proportionalLabel(rep(lseq[i], 2), c(1.2, 1.3), text = FALSE, adj = c(0.5, 0.5), xpd = NA, type = 'l', col = rev(eqCols)[i], lwd = 2)
    }
    LoLinR::proportionalLabel(c(0.24, 0.76, 0.76, 0.24, 0.24), c(1.19, 1.19, 1.31, 1.31, 1.19), text = FALSE, adj = c(0.5, 0.5), xpd = NA, type = 'l')
    LoLinR::proportionalLabel(0.5, 1.1, expression(paste('Transfer efficiency, ' * epsilon, sep = '')), adj = c(0.5, 0.5), xpd = NA)
    
    allEffsCut  <-  allEffs[allEffs <= tippEff[2] + (tippEff[2] - min(allEffs))]
    allRgbs     <-  RGBcomponents(eqCols)
    lseq        <-  c(0.25, 0.37, 0.5, 0.62, 0.75)
    tpos        <-  round(seq(1, length(eqCols), length.out = length(lseq)))
    for (k in seq_along(lseq)) {
        rgbWanted  <-  RGBcomponents(cols[tpos[k]])
        effPos     <-  which.min(apply(rgbWanted - allRgbs, 2, rgbDistance))
        LoLinR::proportionalLabel(rep(lseq[6 - k], 2), c(1.31, 1.33), text = FALSE, adj = c(0.5, 0.5), xpd = NA, type = 'l')
        LoLinR::proportionalLabel(lseq[6 - k], 1.34, round(allEffsCut[effPos], 2), adj = c(0.5, 0), xpd = NA, cex = 0.7)
    }
    
    lpos     <-  c(0.4, 0.3, 0.2)
    for (j in seq_along(tippEff)) {
        LoLinR::proportionalLabel(c(0.74, 0.81), rep(lpos[j], 2), text = FALSE, adj = c(0, 0.5), xpd = NA, type = 'l', lty = ltys[j])
        LoLinR::proportionalLabel(0.83, lpos[j], substitute(italic(epsilon) == b, list(b = unique(round(tippEff[j], 2)))), adj = c(0, 0.5), xpd = NA, cex = 1)
    }
}

##########
# APPENDIX
##########
makeFigureS1  <-  function (dest, ...) {
    toPdf(figS1(...), dest, width = 6, height = 6)
    extrafont::embed_fonts(dest)
}

figS1    <-  function (growthFishBase) {
    mOptMStarRatio  <-  function (b) {
        (1 - 1 / b)^b
    }
    vec  <-  tapply(growthFishBase$b, growthFishBase$species, unique)
    plot(NA, xlab = substitute('Length-mass'~italic('b')~'parameter'), ylab = 'Frequency', xlim = c(2.5, 3.5), ylim = c(0, 120), axes = FALSE)
    plotPreset()
    hist(vec, col = 'black', border = 'white', las = 1, breaks = seq(range(vec)[1], range(vec)[2], length.out = 15), main = '', xlim = c(2.5, 3.5), ylim = c(0, 100), add = TRUE, axes = FALSE, bg = NULL)
    axis(1)
    axis(2, las = 1)
    lines(rep(median(vec), 2), c(0, 115), lty = 2)
    text(median(vec), 117, format(round(mOptMStarRatio(median(vec)), 2), nsmall = 2), adj = c(0.5, 0.5))
    lines(rep(min(vec), 2), c(0, 20), lty = 2)
    text(min(vec), 22, format(round(mOptMStarRatio(min(vec)), 2), nsmall = 2), adj = c(0.5, 0.5))
    lines(rep(max(vec), 2), c(0, 20), lty = 2)
    text(max(vec), 22, format(round(mOptMStarRatio(max(vec)), 2), nsmall = 2), adj = c(0.5, 0.5))
}

makeFigureS2  <-  function (dest, ...) {
    toPdf(figS2(...), dest, width = 9, height = 10)
    extrafont::embed_fonts(dest)
}

figS2  <-  function (growthFishBase, gfStanOut) {
    par(mfcol = c(9, 7), mar = c(0.3, 1.1, 0.3, 0.1), oma = c(4, 4, 2, 2))
    
    xl        <-  expression(paste('Inverse Temperature, 1/', italic(kT[s]), ' - 1/', italic(kT), ' (eV'^{-1}, ')', sep = ''))
    yl        <-  expression(paste('ln ', italic(G['opt']), ' @ 1 g (g d'^{-1}, ')'), sep = '')
    families  <-  sort(unique(growthFishBase$family))
    n         <-  0 #start counting  
    
    for (i in families) {
        n    <-  n + 1
        dat  <-  growthFishBase[growthFishBase$family == i, ]
        
        if (n %in% 1:9) {
            yaxt  <-  's'
        } else {
            yaxt  <-  'n'
        }
        if (n %in% c(9, 18, 27, 36, 45, 54)){
            xaxt  <-  's'
        } else {
            xaxt  <-  'n'
        }
        slope     <-  gfStanOut['lnMopt', 'Estimate'] + gfStanOut[paste0('lnMopt_', unique(dat$family)), 'Estimate']
        response  <-  dat$lnMaxGrowths - slope * dat$lnMopt
        plot(response ~ dat$invKT, xlab = '', ylab = '', las = 1, pch = 21, col = 'grey80', bg = LoLinR::transparentColor('grey80', .5), xaxt = xaxt, yaxt = yaxt, ylim = c(-8, 2), xlim = c(-3, 3))
        
        Int    <-  gfStanOut['Intercept', 'Estimate'] + gfStanOut[paste0('Intercept_', unique(dat$family)), 'Estimate']
        slope  <-  gfStanOut['invKT', 'Estimate'] + gfStanOut[paste0('invKT_', unique(dat$family)), 'Estimate']
        
        lines(range(dat$invKT), Int + slope * range(dat$invKT), lty = 2, col = 'black', lwd = 1)
        LoLinR::proportionalLabel(0.05, 0.85, unique(dat$family), cex = 1, adj = c(0, 0.5))
    }
    mtext(xl, side = 1, line = 2.8, outer = TRUE, cex = 1.2)
    mtext(yl, side = 2, line = 1.7, outer = TRUE, cex = 1.2)
}
