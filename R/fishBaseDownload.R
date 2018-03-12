###################################################################
# THIS IS THE CODE USED TO DOWNLOAD THE FISHBASE DATASET
# IT CAN BE USED TO UPDATE THE DATA IN CASE THERE IS A NEW VERSION
# OF THE DATASETS AVAILABLE ON FISHBASE, PARTICULARLY BECAUSE SOME
# LINKS MAY BECOME DEPRECATED OVER TIME. NOTICE THOUGH THAT BY
# DOWNLOADING UPDATED VERSIONS OF THESE DATASETS FROM FISHBASE,
# THE FINAL ANALYSES MAY LIKELY YIELD DIFFERENT RESULTS BECAUSE
# THE INPUT DATA WILL HAVE CHANGED. ALSO IT IS POSSIBLE THAT THE
# FISHBASE SYSTEM WILL GO DOWN FOR MAINTENANCE WHILE THE DOWNLOAD
# IS HAPPENING, SO KEEP AN EYE OUT BECAUSE YOU MAY HAVE TO RUN THIS
# A FEW TIMES TO GET IT TO WORK. FASTER INTERNET CONNECTIONS HELP.
# PLEASE CONTACT THE AUTHOR VIA THE ISSUES PAGE ON 
# https://github.com/dbarneche/fishgrowth/issues
# SHOULD YOU HAVE PROBLEMS USING THIS CODE.
###################################################################
downloadFishBaseGrowthData  <-  function (fishBaseHtmlFile, verbose = TRUE, ...) {
    taxonInfo   <-  giveTaxonTable(fishBaseHtmlFile)
    gtab  <-  data.frame()
    for (i in 1:nrow(taxonInfo)) {
        gtab  <-  rbind(gtab, getGrowthData(taxonInfo[i, ], verbose = verbose, ...))
    }
    gtab
}

downloadFishBaseLwData  <-  function (fishBaseHtmlFile, verbose = TRUE) {
    taxonInfo  <-  giveTaxonTable(fishBaseHtmlFile)
    lwData  <-  data.frame()
    for (i in 1:nrow(taxonInfo)) {
        lwData  <-  rbind(lwData, getLwData(taxonInfo[i, ], verbose = verbose))
    }
    cleanLwData(lwData)
}

giveTaxonTable  <-  function (pathToHtmlTable) {
    pagetree  <-  XML::htmlTreeParse(pathToHtmlTable, error = function (...){})
    urltable  <-  XML::xmlToList(pagetree$children$html$children$body$children$table$children[[2]])
    
    if (class(urltable) == 'list') {
        data.frame(species  =  unname(sapply(urltable, function (x) unlist(x[1])[1])),
                   urls     =  unname(sapply(urltable, function (x) unlist(x[1])[2])),
                   family   =  unname(sapply(urltable, function (x) unlist(x[3])[1])),
                   stringsAsFactors = FALSE)
    } else if (class(urltable) == 'matrix') {
        data.frame(species  =  unname(apply(urltable, 2, function (x) unlist(x[1])[1])),
                   urls     =  unname(apply(urltable, 2, function (x) unlist(x[1])[2])),
                   family   =  unname(apply(urltable, 2, function (x) unlist(x[3])[1])),
                   stringsAsFactors = FALSE)
    }
}

getGrowthData  <-  function (data, verbose = TRUE, ...) {
    if (verbose) {
        cat('Extracting growth data for:', data$species, '\n')
    }
    sptable         <-  consistentTree(data$urls, verbose = verbose)
    htmlTab         <-  sptable %>% rvest::html_nodes('table')
    while (length(htmlTab) == 0) {
        message('Empty page: re-attempting to download length-weight data for:', data$species, '\n')
        sptable  <-  consistentTree(data$urls, verbose = verbose)
        htmlTab  <-  sptable %>% rvest::html_nodes('table')
    }
    tabPos          <-  grep('Questionable', htmlTab)
    filteredTab     <-  htmlTab %>% magrittr::extract2(tabPos) 
    gtab            <-  (filteredTab %>% rvest::html_table(header = TRUE))[, -1]
    names(gtab)     <-  c('length_inf_cm', 'length_type', 'k_1_over_y', 't_not_years', 'sex', 'm_1_over_y', 'temperature', 'l_m', 'theta', 'country', 'locality', 'questionable', 'captive')
    gtab$species    <-  data$species
    gtab$family     <-  data$family
    gtab$length_inf_cm   <-  as.numeric(gsub(',', '', gtab$length_inf_cm)
    dataTypeHtml    <-  addFishBasePrefix(filteredTab %>% html_nodes('a') %>% html_attr('href'), ...)
    gtab$data_type  <-  sapply(dataTypeHtml, extractDataType, species = data$species, verbose = verbose)
    gtab
}

consistentTree  <-  function (...) {
    sptable  <-  readAndParseTree(...)
    while (inherits(sptable, 'try-error')) {
        sptable  <-  readAndParseTree(...)
    }
    sptable
}

readAndParseTree  <-  function (url, ...) {
    n  <-  1
    cat('Attempt #', n, ' to open: ', url, '\n\n', sep = '')
    sptree  <-  openReadCloseConnection(url, ...)
    while (inherits(sptree, 'try-error')) {
        n  <-  n + 1
        cat('Attempt #', n, ' to open: ', url, '\n\n', sep = '')
        sptree  <-  openReadCloseConnection(url, ...)
    }
    try(xml2::read_html(sptree, verbose = TRUE), silent = TRUE)
}

openReadCloseConnection  <-  function (url, verbose) {
    try(RCurl::getURL(url = url, followLocation = TRUE, .opts = list(timeout = 20, maxredirs = 2, verbose = verbose)), silent = TRUE)
}

addFishBasePrefix  <-  function (sufix, mirror = 'de') {
    sprintf('http://www.fishbase.%s%s', mirror, sufix)
}

extractDataType  <-  function (dataTypeHtml, species, verbose = TRUE) {
    if (verbose) {
        cat('Extracting growth data type for:', species, '\n')
    }
    dataType       <-  consistentTree(dataTypeHtml, verbose = verbose)
    htmlTab        <-  dataType %>% rvest::html_nodes('table')
    while (length(htmlTab) == 0) {
        message('Empty page: re-attempting to download length-weight data for:', species, '\n')
        sptable  <-  consistentTree(dataTypeHtml, verbose = verbose)
        htmlTab  <-  sptable %>% rvest::html_nodes('table')
    }
    tabPos         <-  grep('Data Type', htmlTab)
    filteredTab    <-  htmlTab %>% magrittr::extract2(tabPos) %>% rvest::html_children()
    tabPos         <-  grep('Data Type', filteredTab)
    (filteredTab %>% magrittr::extract2(tabPos) %>% rvest::html_children() %>% rvest::html_text(trim = TRUE))[2]
}

getLwData  <-  function (data, verbose = TRUE) {
    if (verbose) {
        cat('Extracting length-weight data for:', data$species, '\n')
    }
    sptable       <-  consistentTree(data$urls, verbose = verbose)
    htmlTab       <-  sptable %>% rvest::html_nodes('table')
    while (length(htmlTab) == 0) {
        message('Empty page: re-attempting to download length-weight data for:', data$species, '\n')
        sptable  <-  consistentTree(data$urls, verbose = verbose)
        htmlTab  <-  sptable %>% rvest::html_nodes('table')
    }
    tabPos        <-  grep('Score', htmlTab)
    filteredTab   <-  htmlTab %>% magrittr::extract2(tabPos)
    scoresName    <-  filteredTab %>% rvest::html_nodes('input') %>% rvest::html_attr(name = 'name')
    scoresVec     <-  filteredTab %>% rvest::html_nodes('input') %>% rvest::html_attr(name = 'value')
    scoresVec     <-  scoresVec[scoresName == 'score[]']
    gtab          <-  filteredTab %>% rvest::html_table(header = TRUE)
    names(gtab)   <-  c('score', 'a', 'b', 'doubtful', 'sex', 'length_cm', 'length_type', 'r2', 'sd_b', 'sd_log10_a', 'n', 'country', 'locality')
    gtab$score    <-  as.numeric(scoresVec)
    gtab$species  <-  data$species
    gtab$family   <-  data$family
    gtab
}

cleanLwData  <-  function (lwData) {
    for (j in c('score', 'a', 'b', 'sd_b', 'sd_log10_a', 'n')) {
        lwData[[j]]  <-  as.numeric(lwData[[j]])
    }
    
    re               <-  '([[:alpha:].]+)[[:space:]]'
    lwData$doubtful  <-  tolower(gsub(re, '\\1', lwData$doubtful))
    lwData$doubtful  <-  gsub('[[:space:]]', NA, lwData$doubtful)
    
    lwData$sex       <-  tolower(lwData$sex)
    
    lwData$length_cm      <-  fixASCII(lwData$length_cm)
    rangeList             <-  sapply(lwData$length_cm, strsplit, split = '[-][[:space:].]')
    lwData$length_cm_min  <-  sapply(rangeList, evaluateAndReturnRange, pos = 1)
    lwData$length_cm_max  <-  sapply(rangeList, evaluateAndReturnRange, pos = 2)
    
    re                  <-  '([[:alpha:].]+)[[:space:]]'
    lwData$length_type  <-  toupper(gsub(re, '\\1', lwData$length_type))
    lwData$length_type[lwData$length_type == 'NA']  <-  NA
    lwData$length_type  <-  gsub('[[:space:]]', NA, lwData$length_type)
    
    lwData$r2  <-  as.numeric(gsub('&nbsp', '', lwData$r2))
    
    lwData
}

fixASCII  <-  function (charVector) {
    Encoding(charVector)  <-  'latin1'
    iconv(charVector, 'latin1', 'ASCII', sub = '')
}

evaluateAndReturnRange  <-  function (x, pos) {
    lenX  <-  length(x) == 2
    if (lenX) {
        x  <-  gsub(',', '', x)
        as.numeric(x[pos])
    } else {
        NA
    }
}

getRateTraits  <-  function (fishBaseHtmlFile, verbose = TRUE, ...) {
    taxonInfo   <-  giveTaxonTable(fishBaseHtmlFile)
    dietURLs    <-  returnDietURLs(taxonInfo$species, ...)
    dat         <-  data.frame(species = taxonInfo$species,
                               habitat = NA,
                               stringsAsFactors = FALSE)
    for (k in seq_along(dietURLs)) {
        if (verbose) {
            cat('Extracting habitat data for:', dat$species[k], '\n')
        }
        sptable         <-  consistentTree(url = dietURLs[[k]], verbose = verbose)
        dat$habitat[k]  <-  extractAndCleanHabitatInfo(sptable)
    }
    dat
}

returnDietURLs  <-  function (species, mirror = 'de') {
    dietLinks  <-  sapply(species, prepareDietURL)
    sprintf('http://fishbase.%s/summary/%s.html', mirror, dietLinks)
}

prepareDietURL  <-  function (species) {
    speciesPreLink  <-  unlist(strsplit(species, ' '))
    if (length(speciesPreLink) == 2) {
        paste0(speciesPreLink[1], '-', speciesPreLink[2])
    } else {
        paste0(speciesPreLink[1], '-', speciesPreLink[2], '+', speciesPreLink[3])
    }
}

extractAndCleanHabitatInfo  <-  function (...) {
    habitatInfo  <-  extractHabitatInfo(...)
    unname(cleanHabitatInfo(habitatInfo))
}

extractHabitatInfo  <-  function (sptable) {
    blocks        <-  sptable %>% rvest::html_nodes('div')
    blcPos        <-  grep('ss-container', blocks)[1]
    mainBlock     <-  blocks %>% magrittr::extract2(blcPos) %>% rvest::html_nodes('div')
    envPos        <-  grep('ss-main', mainBlock)[1]
    envBlock      <-  mainBlock %>% magrittr::extract2(envPos) %>% rvest::html_children()
    habitatPos    <-  grep('Environment / Climate / Range', envBlock)[1]
    envBlock %>% magrittr::extract(habitatPos+1) %>% rvest::html_children() %>% rvest::html_text(trim = TRUE)
}

cleanHabitatInfo  <-  function (habitatInfo) {
    if (is.na(habitatInfo)) {
        NA
    } else {
        cleanHabitat  <-  sub(' (Ref.', '', habitatInfo, fixed = TRUE)
        standardiseHabitatInfo(cleanHabitat)
    }
}

standardiseHabitatInfo  <-  function (habitatInfo) {
    vec  <-  tolower(gsub('[[:space:].]', '', strsplit(habitatInfo, ';')[[1]]))
    paste0(vec[vec == 'marine' | vec == 'brackish' | vec == 'freshwater' | vec == 'reef-associated' | vec == 'benthopelagic' | vec == 'oceanodromous'], collapse = '; ')
}

prepareFishBaseUrl  <-  function (species) {
    paste0('http://www.fishbase.de/summary/', sub(' ', '-', species), '.html')
}

getTrophicLevelData  <-  function (species, verbose = TRUE) {
    if (species == 'Abramis barellus') {
        sppUrl   <-  prepareFishBaseUrl('Ballerus ballerus')
    } else {
        sppUrl   <-  prepareFishBaseUrl(species)
    }
    sptable  <-  consistentTree(sppUrl, verbose = verbose)
    htmlTab  <-  sptable %>% rvest::html_nodes('div')
    tabPos   <-  grep('Trophic Level', htmlTab, fixed = TRUE)
    if (length(tabPos) == 0) {
        message(species, ' does not contain Trophic Level, consider adding it manually\n')
        data  <-  data.frame(species = species, trophicLevel = NA, se = NA, stringsAsFactors = FALSE)
    } else {
        tlInfo   <-  fixASCII(rvest::html_text(htmlTab[tabPos[length(tabPos)]], trim = TRUE))
        obs      <-  strsplit(tlInfo, '):', fixed = TRUE)[[1]][2]
        numbers  <-  as.numeric(unlist(regmatches(obs, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", obs))))
        data     <-  data.frame(species = species, trophicLevel = numbers[1], se = numbers[2], stringsAsFactors = FALSE)
    }
    data
}

getAspectRatioData  <-  function (species, verbose = TRUE) {
    if (species == 'Abramis barellus') {
        sppUrl   <-  prepareFishBaseUrl('Ballerus ballerus')
    } else {
        sppUrl   <-  prepareFishBaseUrl(species)
    }
    sptable  <-  consistentTree(sppUrl, verbose = verbose)
    htmlTab  <-  sptable %>% rvest::html_nodes('div')
    tag      <-  'physiology/MorphMetList.php?ID='
    tabPos   <-  grep(tag, htmlTab, fixed = TRUE)
    if (length(tabPos) == 0) {
        message(species, ' does not contain Aspect ratio, consider adding it manually\n')
        data  <-  data.frame(species = species, n = NA, aspectRatio_mean = NA, aspectRatio_sd = NA, aspectRatio_min = NA, aspectRatio_max = NA, stringsAsFactors = FALSE)
    } else {
        arTabLink   <-  htmlTab[tabPos[length(tabPos)]] %>% rvest::html_nodes('a') %>% rvest::html_attr('href') %>% grep(pattern = tag, fixed = TRUE, value = TRUE) %>% sub(pattern = '../', replacement = '/', fixed = TRUE) %>% addFishBasePrefix()
        arTab    <-  consistentTree(arTabLink, verbose = verbose) %>% rvest::html_table()
        if (nrow(arTab[[1]]) > 0) {
            arVal    <-  arTab[[1]][['Aspect ratio']]
            data     <-  data.frame(species = species, n = length(arVal), aspectRatio_mean = mean(arVal), aspectRatio_sd = sd(arVal), aspectRatio_min = min(arVal), aspectRatio_max = max(arVal), stringsAsFactors = FALSE)            
        } else {
            data  <-  data.frame(species = species, n = NA, aspectRatio_mean = NA, aspectRatio_sd = NA, aspectRatio_min = NA, aspectRatio_max = NA, stringsAsFactors = FALSE)
        }
    }
    data
}

library(XML)
library(xml2)
library(RCurl)
library(rvest)
library(dplyr)
library(magrittr)

# the code below was run using the following software versions:
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] magrittr_1.5   dplyr_0.7.1    rvest_0.3.2    RCurl_1.95-4.8 bitops_1.0-6   xml2_1.1.1     XML_3.98-1.9  

# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.11     assertthat_0.2.0 R6_2.2.2         httr_1.2.1       rlang_0.1.1      bindrcpp_0.2     glue_1.1.1       compiler_3.4.1   pkgconfig_2.0.1  bindr_0.1        tibble_1.3.3    

fishBaseGrowthData   <-  downloadFishBaseGrowthData('data/fishBaseData/growth_parameters.html')
write.csv(fishBaseGrowthData, 'data/fishBaseData/fishBaseGrowthData.csv', row.names = FALSE)

fishBaseLwData       <-  downloadFishBaseLwData('data/fishBaseData/lw_parameters.html')
write.csv(fishBaseLwData, 'data/fishBaseData/fishBaseLwData.csv', row.names = FALSE)

fishBaseHabitatData  <-  getRateTraits('data/fishBaseData/growth_parameters.html', verbose = TRUE)
write.csv(fishBaseHabitatData, 'data/fishBaseData/fishBaseHabitatData.csv', row.names = FALSE)

# trophic level and aspect ratio
growthFishBase    <-  remake::fetch('growthFishBase')
growthEmData      <-  remake::fetch('growthEmData')
metabolicRates    <-  remake::fetch('metabolicRates')
speciesList       <-  unique(c(growthFishBase$species, growthEmData$species, metabolicRates$species))
trophicLevelData  <-  plyr::ldply(speciesList, getTrophicLevelData)
trophicLevelData[trophicLevelData$trophicLevel > 5  & !is.na(trophicLevelData$trophicLevel), 2:3]  <-  NA # no data
# note that data from species containing NA were manually inserted from closely related species -- this is noted in the obs field
write.csv(trophicLevelData, 'data/fishBaseData/trophicLevelData.csv', row.names = FALSE, quote = FALSE)

aspectRatioData           <-  plyr::ldply(speciesList, getAspectRatioData)
# note that data from species containing NA were manually inserted from closely related species -- this is noted in the obs field
write.csv(aspectRatioData, 'data/fishBaseData/aspectRatioData.csv', row.names = FALSE, quote = FALSE)
