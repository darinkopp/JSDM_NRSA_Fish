# HMSC model fitting script 
# this is used to fit all models

rm(list=ls())
gc()

selectFish <- function(All_NRSA_Fish = fishcnt, 
                       All_NRSA_Sites = siteinfo,
                       VPU = "12", 
                       minOccurrence = 5, 
                       Extend = F){
  
  siteList <- All_NRSA_Sites %>%
    filter(vpu_use==VPU) %>% 
    data.frame(row.names = .$UID) %>%
    dplyr::select(c("SITE_ID", "VISIT_NO", "YEAR", "UID"))
  
  # count occurrences in region
  InRegion.Spp <- All_NRSA_Fish %>%
    filter(UID %in% siteList$UID) %>%
    group_by(FName) %>%
    count(name = "UID") %>%
    filter(UID >= minOccurrence)
  
  fishList <- All_NRSA_Fish %>% 
    filter(UID %in% siteList$UID & FName %in% InRegion.Spp$FName)
  
  # list of species and sites within a network
  InRegion.uid <- unique(fishList$UID)
  
  # The average assemblage richness of sites within the region
  # for example if average richness is 4 then 4 species tend 
  # to co-occur
  MeanR <- fishList %>%
    group_by(UID)%>%
    summarize(rich = n_distinct(FName))%>%
    summarize(round(mean(rich))) %>%
    as.numeric()
  
  # Sites from other regions that have similar assemblages to focal region. 
  # Specifically a site outside th focal region shares a minimum of species 
  # equal to the mean richness of the sites within the focal region.
  addSites <- All_NRSA_Fish %>%
    filter(FName %in% InRegion.Spp$FName & !UID %in% InRegion.uid) %>%
    dplyr::count(UID) %>%
    mutate(Keep = n > MeanR) %>%
    filter(Keep)
  
  addFish <- All_NRSA_Fish %>%
    filter(UID %in% addSites$UID)
  
  addSites <- addFish %>%
    dplyr::select(c(UID, FName, TOTAL)) %>%
    filter(FName%in%InRegion.Spp$FName) %>%
    pivot_wider(names_from = FName, 
                values_from = TOTAL,
                values_fn = ~ ifelse(sum(.x) > 0, 1, 0),
                values_fill = 0)
  
  fishList <- fishList %>%
    dplyr::select(c(UID, FName, TOTAL)) %>%
    filter(FName%in%InRegion.Spp$FName) %>%
    pivot_wider(names_from = FName, 
                values_from = TOTAL,
                values_fn = ~ ifelse(sum(.x) > 0, 1, 0),
                values_fill = 0)
  if(Extend){
    x<-list(InRegion = fishList, OutRegion = addSites) 
  } else {
    x <- fishList
  }
  return(x)
}

#apply logit transformation to place traits simplex
#on normal scale
logit <- function(p) {
  return(log(p / (1 - p)))
}

# load packages
library(Hmsc)
library(ape)
library(tidyverse)

parent.dir <- "C:/Users/DKOPP/OneDrive - Environmental Protection Agency (EPA)/HMSC/Manuscript_JSDM_FishModeling/Data_Analysis_Rscripts"
setwd(parent.dir)
#run this once before starting the script


#setwd(paste0(parent.dir,"/",old.wd))
# Outer loop iterate through thinning allows fitted models to be written without 
# waiting on models that need longer to converge unique(siteinfo$vpu_use)


siteinfo <- read.csv("Data/NRSAFish_Sites.csv")
fishcnt <- read.csv("Data/NRSAFish_Counts_AllYears_Complete.csv")
traits <- read.csv("Data/NRSAFish_Traits_Complete.csv")
phylo <- readRDS("Data/NRSAFish_Tree_Complete.rds")

Regions <- c(unique(siteinfo$vpu_use))
timereport <- data.frame()
for (q in c("10A", "10B", "100A", "100B", "1000A", "1000B")){
  #q <- "1A" 
  print(q)
  # model parameters
  thinning <- as.numeric(substr(q, start = 1, (nchar(q)-1)))
  nChains = 6 #number of chains
  # take a total of 1000 samples by recording 5th observation 
  # 1000 samples with thin = 5 is 5000 iterations for each chain
  thin = thinning #number of iterations per sample
  samples = 500 #number of samples 
  transient  = (samples*thin) #50 * thin #additional iterations used for warmup 
  verbose = 50 * thin #print progress
  
  for (roi in Regions){
    # roi <- "01" 
    print(roi)
    
    Taxa <- selectFish(All_NRSA_Fish = fishcnt,
                       All_NRSA_Sites = siteinfo,
                       VPU = roi, 
                       minOccurrence = 10, 
                       Extend = F)
    
    #will have to change if this is a list...
    Y <- Taxa %>%
      column_to_rownames("UID")
    
    XData <- siteinfo%>%
      filter(UID%in%rownames(Y))%>%
      column_to_rownames("UID")%>%
      select(SITE_ID)%>%
      mutate(SITE_ID=as.factor(SITE_ID))
    
    TrData <- traits %>%
      filter(species %in% colnames(Taxa)) %>%
      column_to_rownames("species") %>%
      select(O, P, E)
    
    #add small non-zero value to support logit transformation 
    TrData <- t(apply(TrData, 1, function(x) ifelse(x==0,.Machine$double.eps,x)))
    TrData <- t(apply(TrData, 1, function(x) logit(x/sum(x))))
    TrData <- data.frame(TrData)
    TrData <- TrData[colnames(Y),]
    
    # phylo tree
    tree <- phylo$phylo
    tree$tip.label <- gsub("_", " ", tree$tip.label)
    phyloTree <- drop.tip(tree, setdiff(tree$tip.label, colnames(Taxa)))
    
    # sample design
    studyDesign = data.frame(SAMPLE = factor(rownames(XData)))
    rL = HmscRandomLevel(units = studyDesign$SAMPLE)
    
    # formulas - only including opportunistic scores. If a taxon is not opportunistic 
    # then it is either periodic or equilibrium (both generally favor more stable environments) 
    # opportunist scores also cover the longest range of values 
    TrFormula = ~ O
    
    # specify model 
    model.env <- Hmsc(Y = 1 * (Y > 0),
                      XData = XData,
                      XFormula=~1,
                      TrData = TrData, 
                      TrFormula = TrFormula,
                      phyloTree = phyloTree,
                      studyDesign = studyDesign, 
                      ranLevels = list("SAMPLE" = rL),
                      dist = "probit")
  
    # fit model
    T1 <- Sys.time()
    print(T1)
    model.env <- sampleMcmc(model.env, thin = thin, 
                            samples = samples, 
                            transient = transient,
                            nChains = nChains, 
                            verbose = verbose, 
                            initPar = "fixed effects", #sets initial condition for MCMC using ML (pg 94)
                            nParallel = nChains) # number of cores cannot be greater than the chains   
    T2 <- Sys.time()
    print(T2 - T1)
    
    #save fitted model. 
    saveRDS(model.env, 
            file = paste0(getwd(),
                          "/Data/Model_Data/HMSC_Model_", 
                          roi, "_", q, "_Intercept.rds"))
    #roi="01"
    #q="100B"
    #thin=100
    # # #if A has already run, add it to B to increase sample size
    # if(length(grep("B", q)) > 0){
    #   initRunA <- paste0("Data/Model_Data/HMSC_Model_", roi, "_",
    #                     paste0(thin, "A"), ".rds")
    #   initRunB <- paste0("Data/Model_Data/HMSC_Model_", roi, "_",
    #                     paste0(thin, "B"), ".rds")
    #  
    #   #model.env <- c(model.env, readRDS(initRun))
    # }
    #model.env.sv<-model.env

    # Check Convergence of model, if satisfactory break loop
    # potential scale reduction factor: goal for PSRF < 1.1. This is slightly high
    # but acceptable given model complexity.
    #####
    # model.env <- readRDS(initRunA)
    # model.env <- readRDS(initRunB)
    # model.env <- c(readRDS(initRunA), readRDS(initRunB))
    # model.env<-model.env.sv
    #is <-grep("Intercept",list.files("Data/Model_Data"), value = T)
    #model.env <- readRDS(paste0("Data/Model_Data/",is[3]))
    
    mpost <- convertToCodaObject(model.env, 
                                 spNamesNumbers = c(T, F), 
                                 covNamesNumbers = c(T, F)) 
  
    
    # Calculate convergence diagnostics (psrf)
    psrf.Beta <- gelman.diag(mpost$Beta, multivariate = F)$psrf
    psrf.Gamma <- gelman.diag(mpost$Gamma, multivariate = F)$psrf
    psrf.Rho <- gelman.diag(mpost$Rho, multivariate = F)$psrf
    
    # # Assoc. parms; Note using the random subsample to avoid excessive 
    # # computation (Pg 305) because we want to use these relationships
     tmp <- mpost$Omega[[1]]
     lapply(tmp, dim)
     z <- dim(tmp[[1]])[2]
     if(z > 1000){sel <- sample(1:z, size = 1000)} else {sel <- 1:z}
     for(j in 1:length(tmp)){tmp[[j]] = tmp[[j]][,sel]}
     psrf.Omega <- coda::gelman.diag(tmp, multivariate = F)$psrf
    
    # check diagnostics
    Beta.psrf.UCI = quantile(psrf.Beta[,"Point est."], 0.95) < 1.1
    Gamma.psrf.UCI = quantile(psrf.Gamma[,"Point est."], 0.95) < 1.1
    Rho.psrf.UCI = quantile(psrf.Rho[,"Point est."], 0.95) < 1.1
    Omega.psrf.UCI = quantile(psrf.Omega[,"Point est."], 0.95) < 1.1
  
    print(c(beta = Beta.psrf.UCI, 
            gamma = Gamma.psrf.UCI, 
            rho = Rho.psrf.UCI, 
            omega = Omega.psrf.UCI))
    
    timereport <- rbind(timereport,data.frame(roi, q, time = T2 - T1,
                                              beta = Beta.psrf.UCI,                      
                                              omega = Omega.psrf.UCI))
    
    # break loop if successfully converged
    # update region vector so it is not re-run with more 
    # iterations
    if(all(Beta.psrf.UCI,Gamma.psrf.UCI,Rho.psrf.UCI,Omega.psrf.UCI)){
      
      Regions <- Regions[Regions != roi]
      #break
      }
    }
}



TrData["Aphredoderus sayanus",]
#hist(psrf.Beta)

psrf.Beta[psrf.Beta[,"Point est."] > 1.1,]
#mean(model.env$Y[,"Notropis volucellus"])
a <- lapply(mpost$Beta, function(x) x[,"B[(Intercept), Carpiodes velifer]"])

mean(model.env$Y[,"Ameiurus catus"])
pnorm(-0.8)

coda::traceplot(a, smooth=T)
coda::traceplot(a[2])
lapply(a, min)


plot(density(unlist(poolMcmcChains(a))))