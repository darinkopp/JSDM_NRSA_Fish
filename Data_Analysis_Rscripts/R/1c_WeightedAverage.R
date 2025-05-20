# Calculating weighted averages for fish species
# Using optimos prome package... for weighted average calcualtions

library(parallel)
library(optimos.prime)
library(dplyr)
library(reshape2)


Fish <- read.csv("Data/NRSAFish_Ocurrence.csv")

#add family to occurrence file and change name because of spelling!
# identify family name for interpretation
#FamName <- unique(Fish_Tax[, c("Accepted_Family", "FName")])

# Environmental Data -- env.data
Env <- read.csv("Data/NRSAFish_Environment.csv")
rownames(Env) <- Env$UID

# calculate optimum for several punitively important environmental gradients
grad <- c("COND_RESULT", "DOC_RESULT", "NTL_RESULT", 
          "PTL_RESULT", "PH_RESULT", "PH_RESULT",
          "TMEAN_S", "TURB_RESULT")

# Select Environmental gradients



# identify taxa with at least 10 occurrences at sites with complete data. 
Include <- table(Fish[, "FName"])>10
Include <- names(Include)[Include]
###################################
Fish_PA <- Fish[Fish$FName %in% Include,]

# Calculate optima for each unique species
###########
Fish <- Fish_PA
spp <- unique(Fish$FName)

clust <- makeCluster(detectCores())
clusterExport(clust, c("Fish", "Env", "grad"))
clusterEvalQ(clust, library(optimos.prime))
optimus <- parLapply(cl = clust, spp, function(i){
  
  #Selects presence records for spp i
  #i<-spp[1]
  Pi <- Fish[Fish$FName == i, "UID"]  
  
  # selects and format assembalge at each Pi
  # The function requires taxa x site matrix
  Pi_assem <- Fish[Fish$UID %in% Pi, ]
  Pi_assem <- reshape2::dcast(FName ~ UID,
                              data = Pi_assem, 
                              fill = 0, 
                              value.var = "TOTAL", 
                              fun.aggregate = sum)
  
  #select and format environmental data at Pi
  Pi_env <- Env[colnames(Pi_assem)[-1], c(grad)]
  Pi_env <- data.frame(ENV = grad, t(Pi_env ))
  
  opt <- op_calculate(environmental_df = Pi_env, 
                      species_df = Pi_assem, 
                      isRelAb = FALSE)
  
  return(opt[opt$Species==i,])
  #optimum_calc <- rbind(optimum_calc, opt[opt$Species==i,])
})
stopCluster(clust)
optimus <- do.call(rbind, optimus)
#optimus <- merge(FamName, optimus, by.x = "FName", by.y = "Species")
###################################

names(optimus)
# format output in long format with estimate and High/Low Tolerance
###########
optimus_long <- melt(optimus, id.vars = c("Species"))
optimus_long$variable <- paste(as.character(optimus_long$variable)," ")

optimus_long$CAST <- optimus_long$variable %>%
  strsplit(., " ") %>%
  lapply(., "[[", 2) %>%
  unlist()

optimus_long$CAST[optimus_long$CAST == ""] <- "Optima"

optimus_long$variable <- optimus_long$variable%>%
  strsplit(., " ") %>%
  lapply(., "[[", 1) %>%
  unlist()

optimus_long <- reshape2::dcast(Species+ variable ~ CAST, 
                                data = optimus_long, 
                                value.var = "value")

###################################

windows()
### using weighted-averaging in the contect of SSD
# https://cran.r-project.org/web/packages/ssdtools/vignettes/ssdtools.html

library(ssdtools)
optimus_long[optimus_long$Species == "Oncorhynchus mykiss",]
######
SSD <- optimus_long[optimus_long$variable=="TMEAN_S",]
# see vignette above.

fits_Opt <- ssd_fit_dists(SSD, left = "Optima",
                          dists = c("llogis_llogis", 
                                    "llogis", 
                                    "gamma","lnorm"), 
                          rescale = T)
ssd_gof(fits_Opt)


fits_range <- ssd_fit_dists(SSD, left = "-LOW",
                            right ="-HIGH",
                            dists = c("llogis_llogis"), 
                            rescale = T)

# The value which 95% of taxa have an optima below. This may serve as a 
# threshold because only 5% of taxa have can thrive above this value.
# this might be too high because lots of taxa have optima below this level. 
# if we manage exclusively to this level sever species will be stressed
ssd_hc(fits_Opt, percent = 95, ci = TRUE)

# The value which 5% of taxa have an optima below. This may serve as a 
# threshold because only 5% of taxa may become stressed above this value. 
# this is conservative and assumes that fish with higher conductivity 
# tolerances can survive low salinity conditions. 

ssd_hc(fits_Opt, proportion = 0.05, ci = TRUE)


# do we bracket values: 
# 90% of taxa will have values between 62.1 and 897.0. This doesn't say anything 
# about which taxa would occur naturally or the taxa that may prefer conditions
# within this range, only that beyond this range, we might expect fish 
# assemblages to start struggling
#########################################
windows()
theme_set(theme_bw()) # set plot theme
autoplot(fits_Opt) + 
  ggtitle("Species Sensitivity Distributions for Specific Conductance") +
  scale_colour_ssd()

autoplot(fits_range) + 
  ggtitle("Species Sensitivity Distributions for Specific Conductance") +
  scale_colour_ssd()

# the range of optima for fish species collected by the national rivers and 
# streams assessment
#########
windows()
ggplot(optimus_long, aes(x = Optima))+
  geom_histogram()+
  facet_wrap("variable", scales= "free")

ggplot(optimus_long, aes(x = Optima, y = Accepted_Family))+
  geom_boxplot(color = "black")+
  scale_x_log10()+
  facet_wrap("variable", scales= "free")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 7))
##################################


# Calculate SSR --- "Sensitive species ratio" from Griffiths 2023
######
SSR_Data <- data.frame()
for (i in unique(optimus_long$variable)){
  # i <- "Tmax9120Ws"
  
  # identify quantiles for each environmental variable
  # and group taxa their sensitivity (incremental values 0,1)
  test <- optimus_long[optimus_long$variable == i,]
  quants <- quantile(test$Optima, probs = (seq(0, 1, 0.1)))
  test$SV <- cut(test$Optima, 
                 breaks = quants, labels = (1:10),
                 include.lowest = T)
  
  # this essentially reduces quantiles foen to down to two groups...
  # category below 50% precentile... these are supposed 
  # to be the most sensitive taxa based on weighted averages
  test$SG <- ifelse(as.numeric(test$SV) < 5, "Sen", "Tol")
  
  # Vector for multiplication 
  # negative = sensitive, Positive value = Tolerant
  test$VAL <- setNames(ifelse(test$SG == "Tol", 1, -1), test$FName)
  
  #format occurrence values
  PA <- reshape2::acast(UID ~ FName, 
                        data = Fish_PA, 
                        value.var = "TOTAL", 
                        sum, fill = 0)
  
  # should always be true
  if(!all(colnames(PA) == test$FName)){stop()}
  
  #Multiply PA counts data by sensitivty classes (-1,1)
  # negative values indicate sensitive taxa and positive values
  # are tolerant taxa. 0's are absent;
  SSR <- PA[, test$FName] * test$VAL
  
  # calculate Sensitive Species Ratio using relative abundance
  SSR_calc <- apply(SSR, 1, function(x) sum(abs(x[x < 0]))/sum(abs(x[x!=0])))
  SSR <- data.frame(UID = rownames(PA), SSR_calc)
  
  # merge back with Env Data to check whether ratio changes along gradient
  SSR <- merge(env.data[c("UID", i)], SSR, by = "UID")
  
  SSR <- setNames(SSR, c("UID", "ENV", "SSR_calc"))
  SSR_Data <- rbind(SSR_Data, data.frame(SSR, variable = i))
}

# Griffith grouped SSR values into bins --- not sure why this was done
SSR_Data$SSR_Bin <- split(SSR_Data,SSR_Data$variable) %>%
  lapply(.,function(x) cut(x$SSR_calc, breaks = 10)) %>%
  unlist()
################################

#merge with site data to explore how the sensitive species ratio varies 
#among sites and with condition classes
######
methods <- read.csv("NRSA_FishMethods.csv")
siteList <- read.csv("Site.info.Master_12112023.csv")

siteList <- merge(siteList, methods, by = c("SITE_ID", "VISIT_NO", "YEAR"))
siteList$cycle <- ifelse(siteList$YEAR %in% c(2008,2009), 
                         "0809", ifelse(siteList$YEAR %in% c(2013, 2014), 
                                        "1314", "1819"))

SSR_Data <- merge(siteList, SSR_Data, by = "UID")
###############################

windows()
#plot shows Binned SSR against environmental gradients 
ggplot(SSR_Data, aes(x = ENV, y = SSR_Bin))+
  geom_jitter(alpha = 0.5)+
  geom_boxplot()+
  scale_x_log10()+
  facet_wrap("variable", scales = "free")

#Plot shows the SSR Clac againast RT Master
ggplot(SSR_Data, aes(x = RT_MASTER, y = SSR_calc))+
  geom_jitter(alpha = 0.5)+
  geom_boxplot()+
  #scale_x_log10()+
  facet_wrap("variable", scales = "free")