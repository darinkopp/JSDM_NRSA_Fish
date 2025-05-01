# Harmonizing NRSA Chemical and Physical Habitat Data 
# data downloaded directly from website or NARS IM. 
# Script creates "Data_Raw/NRSA_ChemPhab.csv"

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
# chemical
file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809_WaterChem_SuperWide_alltheNRSA.tab",
          paste0(getwd(),"/Data/NRSA0809_WaterChem_SuperWide_alltheNRSA.tab"))
file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA1314_WaterChem_SuperWide_alltheNRSA.tab",
          paste0(getwd(),"/Data/NRSA1314_WaterChem_SuperWide_alltheNRSA.tab"))
file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA1819_WaterChem_SuperWide_alltheNRSA.tab",
          paste0(getwd(),"/Data/NRSA1819_WaterChem_SuperWide_alltheNRSA.tab"))

# copy physical habitat
file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809-1819_PhabMetrics_alltheNRSA.tab",
          paste0(getwd(),"/Data/NRSA0809-1819_PhabMetrics_alltheNRSA.tab"))

#copy landscape metrics
file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809-1819_landMets_alltheNRSA.tab",
          paste0(getwd(),"/Data/NRSA0809-1819_landMets_alltheNRSA.tab"))
######################################

FishUID <- read.csv("Data/NRSAFish_Counts_AllYears_Complete.csv") %>%
  select(UID) %>%
  distinct()

#Prepare Site data file - Include survey cycle and vpu 
#####
vpus <- nhdplusTools::get_boundaries(type = "vpu")

#siteinfo files 
siteinfo <- paste0("Data/NRSA0809-1819_siteinfo.tab")|>
  read.table(sep = "\t", header = T)|>
  mutate(HAS_FISH = ifelse(UID%in%FishUID$UID,"Y","N"))|>
  mutate(YEAR = substring(DATE_COL, nchar(DATE_COL)-3, nchar(DATE_COL)))|>
  mutate(CYCLE = ifelse(YEAR %in% c(2008, 2009), "0809", 
                        ifelse(YEAR %in% c(2013, 2014), 
                               "1314", "1819")))|> 
  st_as_sf(coords = c("LON_DD83", "LAT_DD83"), crs = 4269) |>
  st_join(vpus, join = st_within)%>%
  mutate(LON_DD83 = st_coordinates(.)[,1],
         LAT_DD83 = st_coordinates(.)[,2])|>
  st_set_geometry(NULL)|>
  mutate(vpu_use = VPUID)


#Combine VPUs
siteinfo$vpu_use[siteinfo$vpu_use == "14"] <- "15"
siteinfo$vpu_use[siteinfo$vpu_use %in% c("03N", "03S", "03W")] <- "03"
siteinfo$vpu_use[siteinfo$vpu_use %in% c("07", "08")] <- "08"
siteinfo$vpu_use[siteinfo$vpu_use %in% c("05", "06")] <- "05"
# this point had funky coordinates but has a comid and fish data (?)
siteinfo[is.na(siteinfo$vpu_use),"vpu_use"]<-
  substring(siteinfo[is.na(siteinfo$vpu_use),"HUC2"],2)
####################################

#write_site data... The has fish column will allow you to used weights
#this table might change because i want to add the data availability 
#lots of Phab Missing? 
write.csv(siteinfo, "Data/NRSAFish_Sites.csv", row.names = F)
################################################################################
# Prepare water quality
################################################################################

#####
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
  filter(UID%in%FishUID$UID)|>
  dplyr::select(all_of(ChemVarNames))

WQ_1314 <-"Data/NRSA1314_WaterChem_SuperWide_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)|>
  filter(UID%in%FishUID$UID)|>
  dplyr::select(all_of(ChemVarNames))

WQ_1819 <- "Data/NRSA1819_WaterChem_SuperWide_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)|>
  filter(UID%in%FishUID$UID)|>
  dplyr::select(all_of(ChemVarNames))
####################################
chem <- do.call(rbind, list(WQ_0809, WQ_1314, WQ_1819))

siteinfo$HAS_CHEM <- ifelse(siteinfo$UID %in% chem$UID,"Y","N")


# check Units
# apply(chem[,grep("UNITS", names(chem), value = T)], 2, unique)

################################################################################        
# Prepare physical habitat
################################################################################

# general Phab data
#####
# LSUB_DMM Log10(Dgm–Geometric Mean Bed Surface Particle Diameter–mm);
# PCT_SA Sand – .06-2 mm (%); 
# PCT_FN Bed Surface % Fines <0.06mm; 
# "LDCBF_G08"	"SUBSTRATE CHARACTERIZATION"	"Log10(Streambed Critical Diameter-at Bankfull -- mm)(PRK 2008)"
# "W1_HAG"	"HUMAN DISTURBANCE"	"Human Agricultural Influence Index(distance-wtd tally of types and presence"	
# "W1_HALL"	"HUMAN DISTURBANCE"	"Human Disturbance Index(distance-wtd tally of types and presence"
# "W1_HNOAG"	"HUMAN DISTURBANCE"	"Human Non-Agricultural Disturbance Index (distance weighted tally of types and presence)"
# "PCT_FAST"	"CHANNEL HABITAT"	"Percent Fast Water Habitat"
# "PCT_POOL"	"CHANNEL HABITAT"	"Pools -- All Types (% of reach)"
# "XWIDTH"	"CHANNEL MORPHOLOGY"	"Mean Wetted Width (m)"	"Mean Wetted Width (m)"
# "XDEPTH_CM"	"INTERPRETATION"	"Mean thalweg depth (cm)
# "XFC_NAT"	"FISH COVER"	"Sum of non-anthropogenic fish areal cover types"
# "XCMGW"	"RIPARIAN VEGETATION"	"Sum of Woody Canopy+Mid+Ground layer areal cover proportion"	

phabVars <- c("UID", "LSUB_DMM", "PCT_SA", "PCT_FN", 
              "LDCBF_G08","W1_HAG","W1_HALL",
              "PCT_FAST", "PCT_POOL","XWIDTH","XDEPTH_CM",
              "XFC_NAT","XCMGW")

PHAB <- "Data/NRSA0809-1819_PhabMetrics_alltheNRSA.tab"|>
  read.table(sep = "\t", header = T)|>
  filter(UID%in%FishUID$UID)|>
  select(phabVars)
  
apply(PHAB, 2, function(s) sum(is.na(s)))

# Nonsensical value for LDCBF, change to NA
phab0809 <- read.csv("Data_Raw/NRSA_0809/phabmed.csv", na.strings="", stringsAsFactors = F) %>% 
  mutate(LDCBF_G08 = ifelse(LDCBF_G08 == "#NAME?", NA, LDCBF_G08))%>%
  dplyr::select(c("UID", "LSUB_DMM", "PCT_SA", "PCT_FN", "LRBS_G08", "LDCBF_G08",
           "W1_HALL", "W1_HAG", "W1_HNOAG", "XCMGW", "XFC_NAT", "XDEPTH_CM",
           "XWIDTH", "RPXDEP_CM", "PCT_POOL", "PCT_FAST"))%>%
  mutate(LDCBF_G08 = as.numeric(LDCBF_G08))

phab1314 <- read.csv("Data_Raw/NRSA_1314/nrsa1314_phabmed_04232019.csv")%>%
  dplyr::select(c("UID", "LSUB_DMM", "PCT_SA", "PCT_FN", "LRBS_G08", "LDCBF_G08",
                "W1_HALL", "W1_HAG", "W1_HNOAG", "XCMGW", "XFC_NAT", "XDEPTH_CM",
                "XWIDTH", "RPXDEP_CM", "PCT_POOL", "PCT_FAST"))

phab1819 <- read.csv("Data_Raw/NRSA_1819/nrsa_1819_physical_habitat_larger_set_of_metrics_-_data.csv")%>%
  dplyr::select(c("UID", "LSUB_DMM", "PCT_SA", "PCT_FN", "LRBS_G08", "LDCBF_G08",
           "W1_HALL", "W1_HAG", "W1_HNOAG", "XCMGW", "XFC_NAT", "XDEPTH_CM",
           "XWIDTH", "RPXDEP_CM", "PCT_POOL", "PCT_FAST"))

phab <- rbind(phab1819, phab1314, phab0809)
#####################################
NRSA_ChemPhab <- Reduce(function(x,y) merge(x,y,by = "UID", all.x = T), 
                        list(siteinfo, chem, phab))

# add "LQLow_cl" data provided by PRK
######
phab1 <- read.csv("Data_Raw/Darin_OE1819.csv", header = T)
phab2 <- read.csv("Data_Raw/Darin_OE8934.csv", header = T)
phab_LQLow <- rbind(phab2[,intersect(names(phab2),names(phab1))],
              phab1[,intersect(names(phab2),names(phab1))])%>%
  dplyr::select(c("SITE_ID", "VISIT_NO", "YEAR", "LQLow_cl")) %>%
  distinct() #removes duplicated record, phab[phab$SITE_ID=="MORM-1002",]
##################################### 
NRSA_ChemPhab <- merge(NRSA_ChemPhab, phab_LQLow, 
                       all.x = T, 
                       by = c("SITE_ID", "VISIT_NO", "YEAR"))

#add calculated habitat complexity 
######
#vars for habitat complexity
pfc.vars <- c("PFC_ALG","PFC_AQM",
              "PFC_BRS","PFC_LWD",
              "PFC_LVT","PFC_OHV",
              "PFC_RCK","PFC_UCB",
              "PFC_HUM") 

pfc <- read.csv("Data_Raw/NRSA2008-2019_Fish_Cover_PhabMetrics.csv")
pfc <- acast(UID ~ PARAMETER, data = pfc, value.var = "RESULT")

# calculations provided by PK - sum of mean proportion of each cover type across 11 transects 
SumAll_PFC = with(data.frame(pfc), PFC_ALG + PFC_AQM + PFC_BRS + PFC_LWD + PFC_LVT + PFC_OHV + PFC_RCK + PFC_UCB + PFC_HUM);                                    
SumAll_noAlgAQM = with(data.frame(pfc), PFC_BRS + PFC_LWD + PFC_LVT + PFC_OHV + PFC_RCK + PFC_UCB + PFC_HUM)                                                 
SumNat_PFC= with(data.frame(pfc), PFC_ALG + PFC_AQM + PFC_BRS + PFC_LWD + PFC_LVT + PFC_OHV + PFC_RCK + PFC_UCB)                                              
SumBig_PFC= with(data.frame(pfc), PFC_LWD + PFC_BRS + PFC_LVT + PFC_RCK + PFC_OHV + PFC_UCB)                                                                  

#create PA
pfc[pfc > 0] <- 1

BetaPFC_ALL = with(data.frame(pfc), PFC_ALG + PFC_AQM + PFC_BRS + PFC_LWD+ PFC_LVT+ PFC_OHV+ PFC_RCK+ PFC_UCB+ PFC_HUM)       
BetaPFC_ALL_noAlgAqm = with(data.frame(pfc),PFC_BRS+ PFC_LWD+ PFC_LVT+ PFC_OHV+ PFC_RCK+ PFC_UCB+ PFC_HUM)                       
BetaPFC_NAT_w_AlgAqm = with(data.frame(pfc),PFC_ALG+ PFC_AQM+ PFC_BRS+ PFC_LWD+ PFC_LVT+ PFC_OHV+ PFC_RCK+ PFC_UCB)          
BetaPFC_NAT_no_Alg =   with(data.frame(pfc),PFC_AQM+ PFC_BRS+ PFC_LWD+ PFC_LVT+ PFC_OHV+ PFC_RCK+ PFC_UCB)                       
BetaPFC_BIG = with(data.frame(pfc),PFC_LWD + PFC_BRS + PFC_LVT + PFC_RCK + PFC_OHV + PFC_UCB)                                         

AlfaXBetaPFC_ALL = SumAll_PFC * BetaPFC_ALL;                                                                                                
AlfaXBetaPFC_NAT = SumNat_PFC * BetaPFC_NAT_w_AlgAqm;                                                                                       

pfc <- data.frame(UID = rownames(pfc),
                  AlfaXBetaPFC_ALL,
                  AlfaXBetaPFC_NAT, 
                  SumAll_PFC,
                  SumAll_noAlgAQM,
                  BetaPFC_ALL,
                  BetaPFC_ALL_noAlgAqm)
##############################
NRSA_ChemPhab <- merge(NRSA_ChemPhab, pfc, by = "UID", all.x=T)


################################################################################
# Prepare GIS
################################################################################

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

################################################################################
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
methods <- read.table("Data_Raw/fishinfo_wide_0809.tab", sep = "\t", header= T)
methods <- methods[,c("SITE_ID", "VISIT_NO", "YEAR", "FISH_PROTOCOL", "SAMPLED_FISH")]
Fish.Method <- rbind(Fish.Method, methods)

# 1314
methods <- read.table("Data_Raw/wide_fishinfo_1314.tab", sep = "\t", header= T)
methods$YEAR <- substring(methods$DATE_COL, nchar(methods$DATE_COL)-3, nchar(methods$DATE_COL))
methods <- methods[,c("SITE_ID", "VISIT_NO", "YEAR", "FISH_PROTOCOL", "SAMPLED_FISH")]
Fish.Method <- rbind(Fish.Method, methods)

# 1819
methods <- read.table("Data_Raw/nrsa1819_fishinfoWide_newUid.tab", sep = "\t", header= T)
methods$YEAR <- substring(methods$DATE_COL, nchar(methods$DATE_COL)-3, nchar(methods$DATE_COL))
methods <- methods[,c("SITE_ID", "VISIT_NO", "YEAR", "FISH_PROTOCOL", "SAMPLED_FISH")]
Fish.Method <- rbind(Fish.Method, methods)

#make protocol names consistent across surveys Wadeable Sites
unique(Fish.Method[, c("FISH_PROTOCOL")]) 
Fish.Method$FISH_PROTOCOL[Fish.Method$FISH_PROTOCOL %in% c("WADE", "SM_NONWADEABLE", "WADEABLE")] <- "SM_WADEABLE"
Fish.Method$FISH_PROTOCOL[Fish.Method$FISH_PROTOCOL %in% c("LGWADE")] <- "LG_WADEABLE"
Fish.Method$FISH_PROTOCOL[Fish.Method$FISH_PROTOCOL %in% c("BOATABLE","BOAT")] <- "LG_NONWADEABLE"
#####################################

NRSA_ChemPhab <- merge(NRSA_ChemPhab, Fish.Method, 
                       by = c("SITE_ID", "VISIT_NO", "YEAR"), 
                       all.x = T)

# Write dataset 
write.csv(NRSA_ChemPhab, "CleanData/NRSA_ChemPhabGIS_TEST.csv", row.names = F)

library(tidyverse)
allTheNRSA <- "O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/"
allTheNRSA_Fishcts <- paste0(allTheNRSA, "NRSA0809-1819_fishCts_alltheNRSA.tab") %>%
  read.table(sep = "\t", header = T)
names(allTheNRSA_Fishcts)
a<-reshape2::acast(UID~TAXA_ID,data=allTheNRSA_Fishcts, value.var = "TOTAL", fill = 0)
TOTAL=data.frame(UID = rownames(a), TOTAL = rowSums(a))

z<-read.csv("CleanData/NRSA_ChemPhabGIS.csv")
z<-merge(z,TOTAL,by="UID")
unique(z$YEAR)
z2<-unique(z[z$YEAR%in%c(2018,2019), c("YEAR","SAMPLED_FISH")])
z2[order(z2$TOTAL),]
table(z[,c("SAMPLED_FISH", "FISH_PROTOCOL")])
