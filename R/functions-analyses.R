####################
# GENERAL STAN SPECS
####################
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###################
# GENERAL FUNCTIONS
###################
readFile  <-  function (path, ...) {
    read.csv(path, header = TRUE, stringsAsFactors = FALSE, ...)
}

lenUn  <-  function (x) {
    length(unique(x))
}

amplitude  <-  function (x) {
    max(x) - min(x)
}

niceThousands  <-  function (number) {
    formatC(number, format = 'd', big.mark = ',')
}

###################
# EXPERIMENTAL DATA
###################
processDataEm  <-  function (data, familyInfo) {
    ################################
    # LOAD AND BREAK INTO 2 DATASETS
    ################################
    data               <-  data[!(data$stage %in% c('adult', 'juvenile')), ]
    growthEmData       <-  data[ ,c('species', 'stage', 'growthRate', 'growthUnit', 'growthMass', 'growthMassUnit', 'growthMassDryWet', 'growthTempKelvin')]
    respirationEmData  <-  data[ ,c('species', 'stage', 'metabolicRate', 'metabolicUnit', 'metabolicMass', 'metabolicMassUnit', 'metabolicMassDryWet', 'metabolicTempKelvin')]
    
    #######################
    # CONVERT MASS TO GRAMS
    #######################
    growthEmData$growthMass_g                                       <-  growthEmData$growthMass
    growthEmData$growthMass_g[growthEmData$growthMassUnit == 'mg']  <-  growthEmData$growthMass[growthEmData$growthMassUnit == 'mg'] / 1e3
    growthEmData$growthMass_g[growthEmData$growthMassUnit == 'ug']  <-  growthEmData$growthMass[growthEmData$growthMassUnit == 'ug'] / 1e6
    
    respirationEmData$metabolicMass_g                                               <-  respirationEmData$metabolicMass
    respirationEmData$metabolicMass_g[respirationEmData$metabolicMassUnit == 'mg']  <-  respirationEmData$metabolicMass[respirationEmData$metabolicMassUnit == 'mg'] / 1e3
    respirationEmData$metabolicMass_g[respirationEmData$metabolicMassUnit == 'ug']  <-  respirationEmData$metabolicMass[respirationEmData$metabolicMassUnit == 'ug'] / 1e6
    
    #############################
    # CONVERT GROWTH TO GRAMS/DAY
    #############################
    growthEmData$growthRate_g_d                                     <-  growthEmData$growthRate
    growthEmData$growthRate_g_d[growthEmData$growthUnit == 'mg_d']  <-  growthEmData$growthRate[growthEmData$growthUnit == 'mg_d'] / 1e3
    growthEmData$growthRate_g_d[growthEmData$growthUnit == 'ug_d']  <-  growthEmData$growthRate[growthEmData$growthUnit == 'ug_d'] / 1e6
    
    ######################################################
    # CONVERT ALL METABOLIC RATES
    # TO mgO2_h, ASSUME O2 density = 1.429 mg O2 / 1 mL O2
    # http://ocean.ices.dk/Tools/UnitConversion.aspx
    ######################################################
    respirationEmData$metabolicRate_mgO2_h                                                 <-  respirationEmData$metabolicRate
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'j_d']       <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'j_d'] / 5 * 1.429 * 1e3 / 24 / 4.184 / 1e3 # (1 L O2 / 5 kcal) * (1.429 g O2 / 1 L O2) * (1000 mg O2 / 1 g O2) * (1 d / 24 h) * (1 kcal / 4.184 kJ) * (1 kJ / 1000 J) = / 5 * 1.429 * 1e3 / 24 / 4.184 / 1e3
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'nlO2_h']    <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'nlO2_h'] / 1e6 * 1.429 # (1 ml O2 / 1e6 nl O2) * (1.429 mg O2 / 1 ml O2) = 1 / 1e6 * 1.429
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'ulO2_h']    <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'ulO2_h'] / 1e3 * 1.429 # (1 ml O2 / 1e3 ul O2) * (1.429 mg O2 / 1 ml O2) = 1 / 1e3 * 1.429
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'mlO2_h']    <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'mlO2_h'] * 1.429 # (1.429 mg O2 / 1 ml O2) = 1 * 1.429
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'ugO2_h']    <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'ugO2_h'] / 1e3 # (1 mg O2 / 1e3 ug O2)
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'nmolO2_h']  <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'nmolO2_h'] / 1e9 * 32 * 1e3  # (1 mol O2 / 1e9 nmol O2) * (32 g O2 / 1 mol O2) * (1e3 mg O2 / 1 g O2) = 1 / 1e9 * 32 * 1e3
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'umolO2_h']  <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'umolO2_h'] / 1e6 * 32 * 1e3 # (1 mol O2 / 1e6 umol O2) * (32 g O2 / 1 mol O2) * (1e3 mg O2 / 1 g O2) = 1 / 1e6 * 32 * 1e3
    respirationEmData$metabolicRate_mgO2_h[respirationEmData$metabolicUnit == 'mgO2_d']    <-  respirationEmData$metabolicRate[respirationEmData$metabolicUnit == 'mgO2_d'] / 24 # (1 d / 24 h)
    
    #####################################
    # CONVERT ALL METABOLIC RATES TO gC_d
    #####################################
    respirationEmData$metabolicRate_gC_d  <-  respirationEmData$metabolicRate_mgO2_h * 0.009 # (24 h / 1 d) * (1 g / 1e3 mg) * (12 g C / 32 g O2) = 24 / 1e3 * 12 / 32 = 0.009
    
    ####################################
    # CONVERT ALL MASSES TO GRAMS OF WET 
    # MASS ASSUMING THAT 1 GROWTH OF WET 
    # MASS HAS 85% OF WATER
    ####################################
    growthEmData$growthMass_g_wet                                                          <-  growthEmData$growthMass_g
    growthEmData$growthMass_g_wet[growthEmData$growthMassDryWet == 'dry']                  <-  growthEmData$growthMass_g[growthEmData$growthMassDryWet == 'dry'] / 0.15
    growthEmData$growthRate_g_d[growthEmData$growthMassDryWet == 'dry']                    <-  growthEmData$growthRate_g_d[growthEmData$growthMassDryWet == 'dry'] / 0.15
    respirationEmData$metabolicMass_g_wet                                                  <-  respirationEmData$metabolicMass_g
    respirationEmData$metabolicMass_g_wet[respirationEmData$metabolicMassDryWet == 'dry']  <-  respirationEmData$metabolicMass_g[respirationEmData$metabolicMassDryWet == 'dry'] / 0.15
    
    ##################################
    # REMOVE DATA WITH NEGATIVE GROWTH
    ##################################
    respirationEmData  <-  respirationEmData[growthEmData$growthRate_g_d > 0, ]
    growthEmData       <-  growthEmData[growthEmData$growthRate_g_d > 0, ]
    
    ################################
    # REMOVE DATA WHERE MASSES DIFER
    ################################
    masses             <-  growthEmData$growthMass_g_wet == respirationEmData$metabolicMass_g_wet
    growthEmData       <-  growthEmData[masses, ]
    respirationEmData  <-  respirationEmData[masses, ]
    
    ###############
    # CALCULATE Em*
    ###############
    # Metabolic rates, B, are given in g C /day, while growth rates, G, are given in g / day, so B/G results g C / g. Em values reported on Moses et al. (2008) AmNat are expressed in Joules / g, so we need to convert g C to Joules. This conversion should be calculated based on the standard enthalpy combustion of glucose, i.e. a biochemical argument which corresponds to metabolism - it includes both the energy to be sequestrated by ATP and the energy lost during this process. (2805 kJ / 1 mole Glucose) * (1 mole of Glucose / 6 moles C) * (1 mole C / 12 g C) = 2805/6/12 approx. 39 kJ per g C
    growthEmData$emStar  <-  respirationEmData$metabolicRate_gC_d * 39e3 / growthEmData$growthRate_g_d
    
    ####################
    # USE Z SCORE FOR Em 
    # TO REMOVE OUTLIERS
    ####################
    capping            <-  abs(scale(log(growthEmData$emStar))[, 1])
    growthEmData       <-  growthEmData[capping <= 2, ] # outside 95% of the distribution
    respirationEmData  <-  respirationEmData[capping <= 2, ] # outside 95% of the distribution
    
    #############################
    # CREATE COLUMNS FOR ANALYSES
    #############################
    growthEmData$sppNum  <-  as.numeric(as.factor(growthEmData$species))
    growthEmData$stageN  <-  as.numeric(as.factor(growthEmData$stage))
    growthEmData$lnRate  <-  log(growthEmData$growthRate_g_d)
    growthEmData$lnMass  <-  log(growthEmData$growthMass_g_wet)
    growthEmData$invKT   <-  1/8.62e-5 * (1 / 288.15 - 1 / growthEmData$growthTempKelvin)
    
    respirationEmData$sppNum  <-  as.numeric(as.factor(respirationEmData$species))
    respirationEmData$stageN  <-  as.numeric(as.factor(respirationEmData$stage))
    respirationEmData$lnRate  <-  log(respirationEmData$metabolicRate_gC_d)
    respirationEmData$lnMass  <-  log(respirationEmData$metabolicMass_g_wet)
    respirationEmData$invKT   <-  1/8.62e-5 * (1 / 288.15 - 1 / respirationEmData$metabolicTempKelvin)
    
    #############################
    # ADD FAMILY INFO TO PAIRED
    # DATA OF GROWTH & respirationEmData
    #############################
    growthEmData$family       <-  familyInfo$family[match(growthEmData$species, familyInfo$species)]
    respirationEmData$family  <-  familyInfo$family[match(respirationEmData$species, familyInfo$species)]
    growthEmData$famNum       <-  as.numeric(as.factor(growthEmData$family))
    respirationEmData$famNum  <-  as.numeric(as.factor(respirationEmData$family))
    growthEmData$tempC        <-  growthEmData$growthTempKelvin - 273.15
    respirationEmData$tempC   <-  respirationEmData$metabolicTempKelvin - 273.15
    
    list(respirationEmData = respirationEmData, growthEmData = growthEmData)
}

extractGrowthEmData  <-  function (emDatabaseList) {
    emDatabaseList$growthEmData
}

extractRespirationEmData  <-  function (emDatabaseList) {
    emDatabaseList$respirationEmData
}

###############
# FISHBASE DATA
###############
mergeFishBaseData  <-  function (fishBaseGrowthData, fishBaseLwData, fishBaseHabitatData) {
    #######################################
    # 1. IMPORT GROWTH AND LENGTH-MASS DATA
    # 2. KEEP a & b VALUES WITH SCORE > 0.9
    # 3. KEEP b VALUES BETWEEN 2.5 & 3.5
    #######################################
    subG   <-  fishBaseGrowthData[, c('species', 'family', 'length_inf_cm', 'k_1_over_y', 'temperature', 'captive', 'locality', 'data_type')]
    subG   <-  subG[complete.cases(subG), ]
    subG   <-  subG[subG$captive == 'No', ] 
    subL   <-  fishBaseLwData[, c('species', 'family', 'score', 'a', 'b', 'length_cm_min', 'length_cm_max')]
    subL   <-  subL[complete.cases(subL), ]
    subL   <-  subL[subL$score > 0.9, ]
    subL   <-  subL[subL$b > 2.5 & subL$b < 3.5, ]
    subG   <-  subG[subG$species %in% subL$species, ]
    subL   <-  subL[subL$species %in% subG$species, ]
    
    ##############################
    # AVERAGE a & b WITHIN SPECIES
    # MERGE WITH GROWTH DATA
    # TRANSFORM GROWTH TO g d^-1
    # INSTEAD OF g yr^-1
    ##############################
    subL$lnA        <-  log(subL$a)
    subL$lnSizeMin  <-  log(subL$length_cm_min)
    subL$lnSizeMax  <-  log(subL$length_cm_max)
    growthFishBase  <-  data.frame(species = sort(unique(subL$species)), lnA = NA, b = NA, lnSizeMin = NA, lnSizeMax = NA, stringsAsFactors = FALSE)
    for (k in 2:5) {
        growthFishBase[[names(growthFishBase)[k]]]  <-  as.numeric(tapply(subL[[names(growthFishBase)[k]]], subL$species, mean))
    }
    growthFishBase               <-  merge(subG, growthFishBase)
    growthFishBase$tempK         <-  growthFishBase$temperature + 273.15
    growthFishBase$length_inf_cm <-  as.numeric(gsub(',', '', growthFishBase$length_inf_cm, fixed = TRUE)) # large species come with comma (e.g. 1,100 cm)
    growthFishBase$lnMaxGrowths  <-  log(gVopt(exp(growthFishBase$lnA), growthFishBase$k_1_over_y, growthFishBase$length_inf_cm, growthFishBase$b) / 365)
    growthFishBase$lnMopt        <-  log(mVopt(exp(growthFishBase$lnA), growthFishBase$length_inf_cm, growthFishBase$b))
    growthFishBase$invKT         <-  1 / 8.62e-5 * (1 / 288.15 - 1 / growthFishBase$tempK)
    growthFishBase$lnLopt        <-  log(growthFishBase$length_inf_cm*(1 - 1 / growthFishBase$b))
    growthFishBase$goodlnMopt    <-  growthFishBase$lnLopt >= growthFishBase$lnSizeMin & growthFishBase$lnLopt <= growthFishBase$lnSizeMax
    growthFishBase               <-  growthFishBase[growthFishBase$goodlnMopt, ]

    ######################################
    # FILTER FOR FAMILIES WITH AT LEAST 5
    # OBSERVATIONS AND 5 DEGREE TEMP RANGE
    ######################################
    fam5                         <-  tapply(growthFishBase$tempK, growthFishBase$family, function (x) c(max(x) - min(x), length(x)))
    fam5                         <-  names(fam5)[sapply(fam5, function (x) x[1] >= 5 & x[2] >= 5)]
    growthFishBase               <-  growthFishBase[growthFishBase$family %in% fam5, ]
    growthFishBase$familyNumber  <-  as.numeric(as.factor(growthFishBase$family))
    
    #####################
    # IMPORT HABITAT DATA
    #####################
    growthFishBase$habitat  <-  tolower(fishBaseHabitatData$habitat[match(growthFishBase$species, fishBaseHabitatData$species)])
    growthFishBase$habitat  <-  sapply(growthFishBase$habitat, function (x) strsplit(x, '; ')[[1]][1])
    
    growthFishBase
}

readAndProcessMetabolicRates  <-  function (path, trophicLevelData, aspectRatioData) {
    metabolicRates             <-  readFile(path)
    metabolicRates             <-  metabolicRates[!is.na(metabolicRates$activity), ]
    metabolicRates             <-  metabolicRates[metabolicRates$stress == 'none specified' | is.na(metabolicRates$stress), ]
    metabolicRates             <-  metabolicRates[!is.na(metabolicRates$rate), c('species', 'family', 'rate', 'weight', 'temp', 'reef', 'activity')]
    metabolicRates             <-  metabolicRates[metabolicRates$weight > 0, ]
    metabolicRates$tempKelvin  <-  metabolicRates$temp + 273.15
    metabolicRates$invKT       <-  1 / 8.62e-5 * (1 / 288.15 - 1 / metabolicRates$tempKelvin)

    # filter for families with at least 5 observations and 5 degree temperature range
    fam5            <-  tapply(metabolicRates$tempKelvin, metabolicRates$family, function (x) c(max(x)-min(x), length(x)))
    fam5            <-  names(fam5)[sapply(fam5, function (x) x[1] >= 5 & x[2] >= 5)]
    metabolicRates  <-  metabolicRates[metabolicRates$family %in% fam5, ]

    # transform rates for g C / m^2 / day
    metabolicRates$lnRate        <-  log(metabolicRates$rate * (metabolicRates$weight / 1000) * 0.009) # (24 h / 1 d) * (1 g / 1e3 mg) * (12 g C / 32 g O2) = 24 / 1e3 * 12 / 32 = 0.009
    metabolicRates$lnMass        <-  log(metabolicRates$weight)
    metabolicRates$trophicLevel  <-  trophicLevelData$trophicLevel[match(metabolicRates$species, trophicLevelData$species)]
    metabolicRates$aspectRatio   <-  aspectRatioData$aspectRatio_mean[match(metabolicRates$species, aspectRatioData$species)]
    metabolicRates
}

subsetMetabolicRates  <-  function (metabolicRates, type) {
    metabolicRates  <-  metabolicRates[metabolicRates$activity == type, ]
    # need to re-filter for families with at least 5 observations and 5 degree temperature range
    fam5            <-  tapply(metabolicRates$tempKelvin, metabolicRates$family, function (x) c(max(x)-min(x), length(x)))
    fam5            <-  names(fam5)[sapply(fam5, function (x) x[1] >= 5 & x[2] >= 5)]
    metabolicRates[metabolicRates$family %in% fam5, ]
}

calculateEmFishBaseData  <-  function (growthFishBase, standardMetRates, brmsOutput, aspectRatioData, trophicLevelData) {
    growthFishBase               <-  growthFishBase[growthFishBase$family %in% standardMetRates$family, ]
    growthFishBase$trophicLevel  <-  trophicLevelData$trophicLevel[match(growthFishBase$species, trophicLevelData$species)]
    growthFishBase$aspectRatio   <-  aspectRatioData$aspectRatio_mean[match(growthFishBase$species, aspectRatioData$species)]
    coefsStdRates                <-  coef(brmsOutput, pars = 'b')$family
    # Metabolic rates, B, are given in g C /day, while growth rates, G, are given in g / day, so B/G results g C / g. Em values reported on Moses et al. (2008) AmNat are expressed in Joules / g, so we need to convert g C to Joules. This conversion should be calculated based on the standard enthalpy combustion of glucose, i.e. a biochemical argument which corresponds to metabolism - it includes both the energy to be sequestrated by ATP and the energy lost during this process. (2805 kJ / 1 mole Glucose) * (1 mole of Glucose / 6 moles C) * (1 mole C / 12 g C) = 2805/6/12 approx. 39 kJ per g C
    alphaI    <-  coefsStdRates[growthFishBase$family, 'Estimate', 'scalingAlpha_Intercept']
    logitErI  <-  coefsStdRates[growthFishBase$family, 'Estimate', 'logitEr_Intercept']
    lnBoTcI   <-  coefsStdRates[growthFishBase$family, 'Estimate', 'lnBoTs_Intercept'] + log(39e3)
    ToptI     <-  coefsStdRates[growthFishBase$family, 'Estimate', 'Topt_Intercept']
    Ei        <-  coefsStdRates[growthFishBase$family, 'Estimate', 'Ei_Intercept']
    trophic   <-  coefsStdRates[growthFishBase$family, 'Estimate', 'trophicLevelSlope_Intercept']

    growthFishBase$lnMetabolicRate_J_d  <-  lnMetabolicRateSchoolfield(lnBoTc = lnBoTcI, lnTrophicLevel = log(growthFishBase$trophicLevel), trophicSlope = trophic, logitEr = logitErI, Ei = Ei, Topt = ToptI, alpha = alphaI, lnMass = growthFishBase$lnMopt, tempKelvin = growthFishBase$tempK)

    growthFishBase$Em_J_g        <-  exp(growthFishBase$lnMetabolicRate_J_d) / exp(growthFishBase$lnMaxGrowths) * goptCorrection(0.77)
    growthFishBase
}

lnMetabolicRateSchoolfield  <-  function (lnBoTc, lnTrophicLevel, trophicSlope, logitEr, Ei, Topt, alpha, lnMass, tempKelvin,  k = 8.62e-5, Tc = 293.15) {
    Er  <-  Ei / (1 + exp(-logitEr))
    lnBoTc + trophicSlope * lnTrophicLevel + alpha * lnMass + Er / k * (1 / Tc - 1 / tempKelvin) - log(1 + exp(Ei / k * (1 / Topt - 1 / tempKelvin)) * (Er / (Ei - Er)))
}

###############
# THEORY MODELS
###############
mMRatio  <-  function (exponent) {
    (1 / exponent)^(1 / (exponent - 1))
}

goptCorrection  <-  function (exponent) {
    (1 - mMRatio(exponent)^(1 - exponent))
}

gVopt  <-  function (a, k, linf, b) {
    a * k * linf * (linf * (1 - 1 / b))^(b - 1)
}

mVopt  <-  function (a, linf, b) {
    a * (linf * (1 - 1 / b))^b
}

######################
# COMPILED GROWTH DATA
######################
instantaneousGrowth  <-  function (m1, m0, t1, t0) {
    (log(m1 / m0) * m1) / (t1 - t0)
}

####################
# STATISTICAL MODELS
####################
runSchoolFieldFitMetabolicRates  <-  function (metabolicRateSubset) {
    # function does not work with numbers, need to create a
    # dummy variable for Boltzmann constant
    metabolicRateSubset$boltzmannK  <-  8.62e-5

    schoolfieldPriors    <-  c(prior(normal(1, 1e6), nlpar = 'lnBoTs'),
                               prior(normal(1, 1e6), nlpar = 'scalingAlpha'),
                               prior(normal(2, 2), nlpar = 'Ei', lb = 1e-10, ub = 5),
                               prior(normal(-1.2, 1), nlpar = 'logitEr'),
                               prior(normal(298, 5), nlpar = 'Topt', lb = 250, ub = 320))

    schoolfieldInits    <-  list(list(scalingAlpha = 0.8, logitEr = -0.2, lnBoTs = -6, Ei = 1, Topt = 295),
                                 list(scalingAlpha = 0.7, logitEr = -1,   lnBoTs = -7, Ei = 2, Topt = 300),
                                 list(scalingAlpha = 0.6, logitEr = -0.7, lnBoTs = -8, Ei = 3, Topt = 305))

    schoolFieldEquation  <-  brms::bf(lnRate ~ lnBoTs + scalingAlpha * lnMass + (Ei / (1 + exp((-1) * logitEr))) * invKT - log(1 + (Ei / (1 + exp((-1) * logitEr))) / (Ei - (Ei / (1 + exp((-1) * logitEr)))) * exp(Ei / boltzmannK * (1 / Topt - 1 / tempKelvin))), lnBoTs ~ 1 + (1|G|family), scalingAlpha ~ 1 + (1|G|family), logitEr ~ 1 + (1|G|family), Topt ~ 1 + (1|G|family), Ei ~ 1, nl = TRUE)

    set.seed(1)
    brms::brm(schoolFieldEquation, data = metabolicRateSubset, family = gaussian(), prior = schoolfieldPriors, sample_prior = TRUE, chains = 3, cores = 3, iter = 3e4, warmup = 1.5e4, inits = schoolfieldInits, control = list(adapt_delta = 0.99, max_treedepth = 15))
}

runSchoolFieldFitStandardRatesFull  <-  function (standardMetRates) {
    # function does not work with numbers, need to create a
    # dummy variable for Boltzmann constant
    standardMetRates$boltzmannK      <-  8.62e-5
    standardMetRates$lnTrophicLevel  <-  log(standardMetRates$trophicLevel)
    standardMetRates$lnAspectRatio1  <-  log(standardMetRates$aspectRatio + 1) # some species are 0

    schoolfieldPriors    <-  c(prior(normal(1, 1e6),  nlpar = 'lnBoTs'),
                               prior(normal(1, 1e6),  nlpar = 'aspectRatioSlope'),
                               prior(normal(1, 1e6),  nlpar = 'trophicLevelSlope'),
                               prior(normal(1, 1e6),  nlpar = 'scalingAlpha'),
                               prior(normal(2, 2),    nlpar = 'Ei', lb = 1e-10, ub = 5),
                               prior(normal(-1.2, 1), nlpar = 'logitEr'),
                               prior(normal(298, 5),  nlpar = 'Topt', lb = 250, ub = 320))

    schoolfieldInits    <-  list(list(aspectRatioSlope = 0.1, trophicLevelSlope = -0.1, scalingAlpha = 0.8, logitEr = -0.2, lnBoTs = -6, Ei = 1, Topt = 295),
                                 list(aspectRatioSlope = 0, trophicLevelSlope = 0.1, scalingAlpha = 0.7, logitEr = -1,   lnBoTs = -7, Ei = 2, Topt = 300),
                                 list(aspectRatioSlope = -0.1, trophicLevelSlope = 0, scalingAlpha = 0.6, logitEr = -0.7, lnBoTs = -8, Ei = 3, Topt = 305))

    schoolFieldEquation  <-  brms::bf(lnRate ~ lnBoTs + aspectRatioSlope * lnAspectRatio1 + trophicLevelSlope * lnTrophicLevel + scalingAlpha * lnMass + (Ei / (1 + exp((-1) * logitEr))) * invKT - log(1 + (Ei / (1 + exp((-1) * logitEr))) / (Ei - (Ei / (1 + exp((-1) * logitEr)))) * exp(Ei / boltzmannK * (1 / Topt - 1 / tempKelvin))), lnBoTs ~ 1 + (1|G|family), scalingAlpha ~ 1 + (1|G|family), logitEr ~ 1 + (1|G|family), Topt ~ 1 + (1|G|family), aspectRatioSlope + trophicLevelSlope + Ei ~ 1, nl = TRUE)

    set.seed(1)
    brms::brm(schoolFieldEquation, data = standardMetRates, family = gaussian(), prior = schoolfieldPriors, sample_prior = TRUE, chains = 3, cores = 3, iter = 3e4, warmup = 1.5e4, inits = schoolfieldInits, control = list(adapt_delta = 0.99, max_treedepth = 15))
}

runSchoolFieldFitStandardRatesFullMinusAspectRatio  <-  function (standardMetRates) {
    # function does not work with numbers, need to create a
    # dummy variable for Boltzmann constant
    standardMetRates$boltzmannK      <-  8.62e-5
    standardMetRates$lnTrophicLevel  <-  log(standardMetRates$trophicLevel)

    schoolfieldPriors    <-  c(prior(normal(1, 1e6),  nlpar = 'lnBoTs'),
                               prior(normal(1, 1e6),  nlpar = 'trophicLevelSlope'),
                               prior(normal(1, 1e6),  nlpar = 'scalingAlpha'),
                               prior(normal(2, 2),    nlpar = 'Ei', lb = 1e-10, ub = 5),
                               prior(normal(-1.2, 1), nlpar = 'logitEr'),
                               prior(normal(298, 5),  nlpar = 'Topt', lb = 250, ub = 320))

    schoolfieldInits    <-  list(list(trophicLevelSlope = -0.1, scalingAlpha = 0.8, logitEr = -0.2, lnBoTs = -6, Ei = 1, Topt = 295),
                                 list(trophicLevelSlope = 0.1, scalingAlpha = 0.7, logitEr = -1,   lnBoTs = -7, Ei = 2, Topt = 300),
                                 list(trophicLevelSlope = 0, scalingAlpha = 0.6, logitEr = -0.7, lnBoTs = -8, Ei = 3, Topt = 305))

    schoolFieldEquation  <-  brms::bf(lnRate ~ lnBoTs + trophicLevelSlope * lnTrophicLevel + scalingAlpha * lnMass + (Ei / (1 + exp((-1) * logitEr))) * invKT - log(1 + (Ei / (1 + exp((-1) * logitEr))) / (Ei - (Ei / (1 + exp((-1) * logitEr)))) * exp(Ei / boltzmannK * (1 / Topt - 1 / tempKelvin))), lnBoTs ~ 1 + (1|G|family), scalingAlpha ~ 1 + (1|G|family), logitEr ~ 1 + (1|G|family), Topt ~ 1 + (1|G|family), trophicLevelSlope + Ei ~ 1, nl = TRUE)

    set.seed(1)
    brms::brm(schoolFieldEquation, data = standardMetRates, family = gaussian(), prior = schoolfieldPriors, sample_prior = TRUE, chains = 3, cores = 3, iter = 3e4, warmup = 1.5e4, inits = schoolfieldInits, control = list(adapt_delta = 0.99, max_treedepth = 15))
}

runModelSelection  <-  function (complexModel, simplerModel) {
    # complex model first, simpler model later
    # _The difference will be positive if the expected predictive accuracy for the second model is higher._ 
    looCompTabIc  <-  brms::LOO(complexModel, simplerModel, pointwise = FALSE, cores = 3)
    looCompTab    <-  loo::compare(looCompTabIc[[1]], looCompTabIc[[2]])
    # because this is a MCMC exercise, we can use the full posterior distribution of 
    # the difference in expected predictive accuracy between both models to calculate a 
    # p-value, i.e. the s.e. distribution approximates a Gaussian because there's a 
    # large number of independent MCMC samples.
    pVal          <-  2 * pnorm(-abs((looCompTab['elpd_diff'] - 0) / looCompTab['se']))
    list(looCompTabIc  =  looCompTabIc,
         looCompTab    =  looCompTab,
         pVal          =  pVal
        )
}

runStanModelsExperimentalData  <-  function (respirationEmData, growthEmData) {
    # function does not work with numbers, need to create a
    # dummy variable for Boltzmann constant
    respirationEmData$boltzmannK   <-  8.62e-5
    growthEmData$boltzmannK        <-  8.62e-5
    respirationEmData$tempKelvin   <-  respirationEmData$metabolicTempKelvin
    growthEmData$tempKelvin        <-  growthEmData$growthTempKelvin

    schoolfieldPriors    <-  c(prior(normal(1, 1e6), nlpar = 'lnBoTs'),
                               prior(normal(1, 1e6), nlpar = 'scalingAlpha'),
                               prior(normal(2, 2), nlpar = 'Ei', lb = 1e-10, ub = 5),
                               prior(normal(-1.2, 1), nlpar = 'logitEr'),
                               prior(normal(298, 5), nlpar = 'Topt', lb = 250, ub = 320))

    schoolfieldInits    <-  list(list(scalingAlpha = 0.8, logitEr = -0.2, lnBoTs = -6, Ei = 1, Topt = 295),
                                 list(scalingAlpha = 0.7, logitEr = -1,   lnBoTs = -7, Ei = 2, Topt = 300),
                                 list(scalingAlpha = 0.6, logitEr = -0.7, lnBoTs = -8, Ei = 3, Topt = 305))

    schoolFieldEquation  <-  brms::bf(lnRate ~ lnBoTs + scalingAlpha * lnMass + (Ei / (1 + exp((-1) * logitEr))) * invKT - log(1 + (Ei / (1 + exp((-1) * logitEr))) / (Ei - (Ei / (1 + exp((-1) * logitEr)))) * exp(Ei / boltzmannK * (1 / Topt - 1 / tempKelvin))), lnBoTs ~ 1 + (1 | species), scalingAlpha + Ei + logitEr + Topt ~ 1, nl = TRUE)
    
    set.seed(1)
    refit      <-  brms::brm(schoolFieldEquation, data = respirationEmData, family = gaussian(), prior = schoolfieldPriors, sample_prior = TRUE, chains = 3, cores = 3, iter = 5e3, warmup = 2.5e3, inits = schoolfieldInits, control = list(adapt_delta = 0.99))
    reStanOut  <-  extractLikeJagsSummaryBrms(refit, groupingVar = 'species', getEr = TRUE)

    set.seed(1)
    gefit      <-  brms::brm(schoolFieldEquation, data = growthEmData, family = gaussian(), prior = schoolfieldPriors, sample_prior = TRUE, chains = 3, cores = 3, iter = 5e3, warmup = 2.5e3, inits = schoolfieldInits, control = list(adapt_delta = 0.99))
    geStanOut  <-  extractLikeJagsSummaryBrms(gefit, groupingVar = 'species', getEr = TRUE)
    
    list(reStanOut = reStanOut, geStanOut = geStanOut)
}

runStanModelsFishBaseData  <-  function (growthFishBase) {
    set.seed(1)
    gffit      <-  brms::brm(lnMaxGrowths ~ lnMopt + invKT + (1 + lnMopt + invKT | family), data = growthFishBase, family = gaussian(), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 5e3, warmup = 2.5e3)
    gfStanOut  <-  extractLikeJagsSummaryBrms(gffit, groupingVar = 'family')
        
    list(gfStanOut = gfStanOut)
}

extractLikeJagsSummaryBrms  <-  function (stanFit, groupingVar, getEr = FALSE) {
    fixed   <-  as.data.frame(summary(stanFit)$fixed)
    random  <-  as.data.frame(summary(stanFit)$random)
    names(random)  <-  names(fixed)
    error   <-  as.data.frame(summary(stanFit)$spec_pars)
    if (getEr) {
        logitErs  <-  brms::posterior_samples(stanFit, pars = 'logitEr_Intercept')[, 1]
        Eis       <-  brms::posterior_samples(stanFit, pars = 'Ei_Intercept')[, 1]
        Ers       <-  Eis / (1 + exp((-1) * logitErs))
        groupRe   <-  ranef(stanFit)[[groupingVar]]
        rownames(groupRe)  <-  paste0(dimnames(groupRe)[[3]], '_', rownames(groupRe))
        colnames(groupRe)[3:4]  <-  c('l-95% CI', 'u-95% CI')
        groupRe   <-  data.frame(groupRe[, , 1],  'Eff.Sample' = NA, 'Rhat' = NA, check.names = FALSE)
        rbind(fixed, data.frame('Estimate' = mean(Ers), 'Est.Error' = sd(Ers) / sqrt(length(Ers)), 'l-95% CI' = quantile(Ers, probs = 0.025), 'u-95% CI' = quantile(Ers, probs = 0.975), 'Eff.Sample' = NA, 'Rhat' = NA, row.names = 'Er', check.names = FALSE), random, error, groupRe)
    } else {
        groupRe  <-  ranef(stanFit)[[groupingVar]]
        ranefs   <-  data.frame()
        for (k in seq_len(dim(groupRe)[3])) {
            dat  <-  groupRe[, , k]
            rownames(dat)  <-  paste0(dimnames(groupRe)[[3]][k], '_', rownames(dat))
            colnames(dat)[3:4]  <-  c('l-95% CI', 'u-95% CI')
            ranefs  <-  rbind(ranefs, dat)
        }
        ranefs   <-  cbind(ranefs,  'Eff.Sample' = NA, 'Rhat' = NA)
        rbind(fixed, random, error, ranefs)
    }
}

extractStanObjects  <-  function (stanOutputList, objName) {
    stanOutputList[[objName]]
}

estimateTempDependenceEm  <-  function (data, response, predictors, randomSlope, randomGroup) {
    eval(parse(text = sprintf('brms::brm(log(%s) ~ %s + (1 + %s | %s), data = data, family = gaussian(), prior = c(prior(normal(1, 2), "b"), prior(normal(3, 3), "Intercept"), prior(student_t(3, 0, 20), "sd"), prior(student_t(3, 0, 20), "sigma")), sample_prior = TRUE, chains = 3, cores = 3, iter = 5e3, warmup = 2.5e3)', response, predictors, randomSlope, randomGroup)))
}

extractFullPosteriorRandom  <-  function (brmsModel, parameterName, groupingLevel) {
    pstSlpFxd  <-  brms::posterior_samples(brmsModel, pars = sprintf('b_%s', parameterName))
    pstSlpRnd  <-  brms::posterior_samples(brmsModel, pars = sprintf(',%s]', parameterName))
    names(pstSlpRnd)  <-  gsub(sprintf(',%s]', parameterName), '', gsub(sprintf('r_%s[', groupingLevel), '', names(pstSlpRnd), fixed = TRUE))

    # get random slope coefficients for caterpillar plot
    postSlopeCoef    <-  data.frame()
    for (k in 1:ncol(pstSlpRnd)) {
        posterior      <-  pstSlpFxd + pstSlpRnd[, k]
        cis95          <-  quantile(posterior[, 1], probs = c(0.025, 0.975))
        postSlopeCoef  <-  rbind(postSlopeCoef, data.frame(group = names(pstSlpRnd)[k], groupNumber = k, mean = mean(posterior[, 1]), lower_95ci = cis95[1], upper_95ci = cis95[2], stringsAsFactors = FALSE, row.names = NULL))
    }
    postSlopeCoef
}
