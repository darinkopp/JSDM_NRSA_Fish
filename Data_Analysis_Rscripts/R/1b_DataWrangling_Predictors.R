# Harmonizing NRSA Chemical and Physical Habitat Data 

rm(list = ls())
gc()

library(reshape2)
library(lubridate)
library(sf)
library(tidyverse)
library(parallel)
library(nhdplusTools)
library(remotes)
library(StreamCatTools)

#copy predictor tables
#############
# # chemical
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809_WaterChem_SuperWide_alltheNRSA.tab",
#           paste0(getwd(),"/Data/NRSA0809_WaterChem_SuperWide_alltheNRSA.tab"))
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA1314_WaterChem_SuperWide_alltheNRSA.tab",
#           paste0(getwd(),"/Data/NRSA1314_WaterChem_SuperWide_alltheNRSA.tab"))
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA1819_WaterChem_SuperWide_alltheNRSA.tab",
#           paste0(getwd(),"/Data/NRSA1819_WaterChem_SuperWide_alltheNRSA.tab"))
# 
# # copy physical habitat
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809-1819_PhabMetrics_alltheNRSA.tab",
#           paste0(getwd(),"/Data/NRSA0809-1819_PhabMetrics_alltheNRSA.tab"))
# 
# #copy landscape metrics
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809-1819_landMets_alltheNRSA.tab",
#           paste0(getwd(),"/Data/NRSA0809-1819_landMets_alltheNRSA.tab"))
######################################

# occurrence 
FishUID <- read.csv("Data/NRSAFish_Ocurrence.csv") %>%
  select(UID) %>%
  distinct()

# siteinfo files 
#########
siteinfo <- paste0("Data/NRSA0809-1819_siteinfo.tab")|>
  read.table(sep = "\t", header = T)|>
  mutate(vpu_use = substring(HUC2,2), 
         YEAR = substring(DATE_COL,7))
  
# Combine VPUs
siteinfo$vpu_use[siteinfo$vpu_use == "14"] <- "15"
siteinfo$vpu_use[siteinfo$vpu_use %in% c("03N", "03S", "03W")] <- "03"
siteinfo$vpu_use[siteinfo$vpu_use %in% c("07", "08")] <- "08"
siteinfo$vpu_use[siteinfo$vpu_use %in% c("05", "06")] <- "05"

siteinfo <- siteinfo%>%
  select("UNIQUE_ID", "UID", 
         "SITE_ID","VISIT_NO",
         "LAT_DD83", "LON_DD83",
         "DATE_COL",
         "DSGN_CYCLE", "AG_ECO9",
         "COMID","HUC2", "WGT_CAT",
         "WGT_TP", "vpu_use", "YEAR")

####################################
anyDuplicated(siteinfo[,c("SITE_ID","VISIT_NO","YEAR")])

###############################################################################
# Methods data obtained from
# L:\Priv\CORFiles\IM-TH007\data\im\nrsaYYYY\raw. From these tables, selected the 
# the protocol used to sample fish "FISH_PROTOCOL", and whether a site was 
# sampled sufficiently "SAMPLED_FISH". This analysis focus on fish collected by 
# electrofishing (majority of the sites) and sampled sufficiently to accurately 
# capture the fish diversity at the site
################################################################################

Fish.Method <- data.frame()
#######
# 0809 
methods <- read.table("Data/fishinfo_wide_0809.tab", sep = "\t", header= T)
methods <- methods[,c("SITE_ID", "VISIT_NO", "YEAR", "FISH_PROTOCOL", "SAMPLED_FISH")]
Fish.Method <- rbind(Fish.Method, methods)

# 1314
methods <- read.table("Data/wide_fishinfo_1314.tab", sep = "\t", header= T)
methods$YEAR <- substring(methods$DATE_COL, nchar(methods$DATE_COL)-3, nchar(methods$DATE_COL))
methods <- methods[,c("SITE_ID", "VISIT_NO", "YEAR", "FISH_PROTOCOL", "SAMPLED_FISH")]
Fish.Method <- rbind(Fish.Method, methods)

# 1819
methods <- read.table("Data/nrsa1819_fishinfoWide_newUid.tab", sep = "\t", header= T)
methods$YEAR <- substring(methods$DATE_COL, nchar(methods$DATE_COL)-3, nchar(methods$DATE_COL))
methods <- methods[,c("SITE_ID", "VISIT_NO", "YEAR", "FISH_PROTOCOL", "SAMPLED_FISH")]
Fish.Method <- rbind(Fish.Method, methods)

#make protocol names consistent across surveys Wadeable Sites
unique(Fish.Method[, c("FISH_PROTOCOL")]) 
Fish.Method$FISH_PROTOCOL[Fish.Method$FISH_PROTOCOL %in% c("WADE", "WADEABLE")] <- "SM_WADEABLE"
Fish.Method$FISH_PROTOCOL[Fish.Method$FISH_PROTOCOL %in% c("SM_NONWADEABLE")] <- "SM_NONWADEABLE"
Fish.Method$FISH_PROTOCOL[Fish.Method$FISH_PROTOCOL %in% c("LGWADE")] <- "LG_WADEABLE"
Fish.Method$FISH_PROTOCOL[Fish.Method$FISH_PROTOCOL %in% c("BOATABLE","BOAT")] <- "LG_NONWADEABLE"
#####################################
anyDuplicated(Fish.Method[,c("SITE_ID", "VISIT_NO", "YEAR")])

# merge siteinfo to the fish methods
EnvData <- merge(siteinfo,Fish.Method,by = c("SITE_ID", "VISIT_NO", "YEAR"), all.x=T)

################################################################################
# Prepare water quality
################################################################################

#########
# units from NRSA site metadata
# Specific Conductance	uS/cm; 
# Sulfate	CHEMW	mg/L; 
# Chloride	CHEMW	mg/L;
# Laboratory measured pH;
# Total Phosphorus	ug/L;
# Total Nitrogen	mg/L;
# DOC	Dissolved Organic Carbon	mg/L

ChemVarNames <- c("UID", 
                  "CHLORIDE_RESULT", "CHLORIDE_UNITS",
                  "SULFATE_RESULT", "SULFATE_UNITS",
                  "COND_RESULT", "COND_UNITS",
                  "DOC_RESULT", "DOC_UNITS",
                  "NTL_RESULT", "NTL_UNITS",
                  "PTL_RESULT", "PTL_UNITS",
                  "PH_RESULT", "PH_UNITS",
                  "TURB_RESULT", "TURB_UNITS")

WQ_0809 <-"Data/NRSA0809_WaterChem_SuperWide_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)|>
  dplyr::select(all_of(ChemVarNames))

WQ_1314 <-"Data/NRSA1314_WaterChem_SuperWide_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)|>
  dplyr::select(all_of(ChemVarNames))

WQ_1819 <- "Data/NRSA1819_WaterChem_SuperWide_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)|>
  dplyr::select(all_of(ChemVarNames))
####################################
chem <- do.call(rbind, list(WQ_0809, WQ_1314, WQ_1819))

# check Units
# apply(chem[,grep("UNITS", names(chem), value = T)], 2, unique)
chem <- chem[,c("UID", grep("RESULT", names(chem), value = T))]

################################################################################        
# Prepare physical habitat
################################################################################

# general fish cover data
#####
# Natural Fish Cover Features
c("PFC_ALG", "PFC_AQM", "PFC_BRS", 
  "PFC_LWD", "PFC_LVT", "PFC_OHV", 
  "PFC_RCK", "PFC_UCB")

# Depth heterogenity
fishHab <- "Data/NRSA0809-1819_PhabMetrics_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)|>
  select(c("UID", 
           "PFC_ALG", "PFC_AQM", 
           "PFC_BRS", "PFC_LWD", 
           "PFC_LVT", "PFC_OHV", 
           "PFC_RCK", "PFC_UCB",
           "XCMGW", "RPVDEP"))
fishHab$FCR <- apply(fishHab[,c(2:9)], 1, function(x) sum(x>0))
  
apply(fishHab, 2, function(s) sum(is.na(s)))

############################
fishHab

################################################################################
# Prepare GIS: Area, Temp, Precip
################################################################################

lsVars <- "Data/NRSA0809-1819_landMets_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)
  
out <- data.frame()
for (i in siteinfo$UID){
  #i<-siteinfo$UID[1]
  
  siy <- siteinfo[siteinfo$UID==i,c("UNIQUE_ID","SITE_ID","YEAR")]
  
  tmp <- lsVars %>%
    filter(UNIQUE_ID == siy$UNIQUE_ID & SITE_ID == siy$SITE_ID)%>%
    mutate(YEAR = siy$YEAR)%>%
    select(c("UNIQUE_ID", "WSAREASQKM",
             paste0("TMEAN_S_", siy$YEAR,"WS"),
             paste0("PSUMPY_", siy$YEAR,"WS"), "YEAR"))%>%
    rename_with(~c("UNIQUE_ID","WSAREASQKM","TMEAN_S","PSUMPY", "YEAR"))
  out <- rbind(out,tmp)         
}
Landscape <- unique(out)
################################
anyDuplicated(Landscape[,c("UNIQUE_ID", "YEAR")])

# merge data sets together
EnvData <- Reduce(function(x,y) merge(x, y, by = "UID", all.x = TRUE), list(EnvData, chem, fishHab))
EnvData <- merge(EnvData,Landscape, by = c("UNIQUE_ID", "YEAR"), all.x =T)

nrow(EnvData)
nrow(siteinfo)

# Write dataset 
write.csv(EnvData, "Data/NRSAFish_Environment.csv", row.names = F)


###############################################################################
# extras to extract additional GIS variables
###############################################################################
library(dataRetrieval)
NRSA_GIS <- NRSA_ChemPhab |>
  dplyr::select(all_of(c("UID", "LAT_DD83", "LON_DD83", 
                         "VPUID", "CYCLE","COMID")))

# Assign streamcat variables initial list
###########
vars <- c("BFI","DamNrmStor", "Elev", "kffact",
          "rddens","popden2010","tmin9120","tmax9120","tmean9120")

# Function limits queries iterate through site groups 
gr <- as.factor(cut(1:nrow(NRSA_GIS),nrow(NRSA_GIS)/100, labels = F))

#sc_get_params("roads")
count <- 1
for(i in vars){
  # i <- "RoadDens"
  sc <- data.frame()  
  for (j in levels(gr)){
    #j <- 1
    sc.tmp <- sc_get_data(metric = i, 
                          aoi = 'catchment, watershed', 
                          comid = NRSA_GIS$COMID[gr==j])
    sc <- rbind(sc, sc.tmp)
  }
  
  #keep catchment area for first variable
  if(count == 1){
    NRSA_GIS <- merge(NRSA_GIS, 
                      unique(sc), 
                      by = "COMID",
                      all.x = T)  
    #drop catchment area 
  } else {
    NRSA_GIS <- merge(NRSA_GIS, 
                      unique(sc[,-c(2,3)]), 
                      by = "COMID",
                      all.x = T)
  }
  
  count<-count + 1
}
######################################
names(NRSA_GIS)

# sv <- NRSA_GIS
# NRSA_GIS <- sv

# Assign NLCD streamcat variables 
###########
nlcd.full <- data.frame()
for (i in unique(NRSA_GIS$CYCLE)){
  #i <- "0809"
  set <- NRSA_GIS[NRSA_GIS$CYCLE==i,c("UID", "COMID")]
  gr <- as.factor(cut(1:nrow(set),
                      nrow(set)/500, 
                      labels = F))
  nlcd <- data.frame()
  for (g in levels(gr)){
    # g <- 1
    
    if(i == "0809"){
      nlcd.tmp <- sc_nlcd(comid = set$COMID[gr==g], year = "2008")  
    }
    
    if(i == "1314"){
      nlcd.tmp <- sc_nlcd(comid = set$COMID[gr==g], year = "2013")
    }
    
    if(i == "1819"){
      nlcd.tmp <- sc_nlcd(comid = set$COMID[gr==g], year = "2019")
    }
    
    nlcd.tmp <- merge(set[gr==g,], nlcd.tmp,by="COMID")
    nlcd <- rbind(nlcd, nlcd.tmp)
  }
  
  #strip Years from names
  names(nlcd) <- gsub("([0-9]+)", "\\2", names(nlcd))
  
  nlcd.full <- rbind(nlcd.full, nlcd)
}

NRSA_GIS <- merge(NRSA_GIS, 
                  subset(nlcd.full, select = -c(COMID, WSAREASQKM, CATAREASQKM)), 
                  by = "UID", 
                  all.x = T)
#################################
names(NRSA_GIS)

# Get slope from NHDPlusV2
########
# NHDplusV2 Data stored locally
nhd.dir <- "L:/Public/dkopp/NHDPLUSV2"
q <- split(NRSA_GIS, NRSA_GIS$VPU)
out <- data.frame()
for (i in 1:length(q)){
  #i = 12
  VPU = unique(q[[i]]$VPU) 
  print(VPU)
  T1 <- Sys.time()
  
  #read NHDPLUS Slope values for temp stations
  #####
  slope <- paste0("NHDPlus", VPU,"/NHDPlusAttributes/elevslope.dbf")
  slope <- grep(slope, list.files(nhd.dir, 
                                  full.names = T, 
                                  recursive=T),
                value = T)
  slope = foreign::read.dbf(slope)
  
  tmp <- slope[slope$COMID %in% q[[i]]$COMID, c("COMID","SLOPE")]
  out <- rbind(out, tmp)
}

NRSA_GIS <- merge(NRSA_GIS, out, by = "COMID", all.x = T) 

NRSA_GIS <- NRSA_GIS |>
  dplyr::select(-c("CYCLE", "VPUID", "LON_DD83", "LAT_DD83")) 
#################################

NRSA_ChemPhab <- merge(NRSA_ChemPhab, 
                       NRSA_GIS, by = "UID", 
                       all.x = T)

#
