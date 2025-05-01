rm(list = ls())
################################################################################
# Organize raw data into usable formats
# Harmonize taxonomic names across Occurrence, Density, traits and 
# phylogeny data sources. The objective here is to create a common name (FName) 
# that can crosswalk across datasets. Usually this is the accepted name 
# from FishBase/AFS but in some instances the synonym that matched the record 
# in phylogeny was used.

# Files produced in CleanData
################################################################################

#install.packages("openxlsx")
#install.packages("U.Taxonstand")
#install.packages("ape")
#devtools::install_github("ecoinfor/U.Taxonstand")
#devtools::install_github("james-thorson/FishLife", dep=T)
#devtools::install_github("jinyizju/U.PhyloMaker")

library(tidyverse)
library(openxlsx) # Read .xls FishBase Names
library(U.Taxonstand) # Taxonomic harmonization 
library("U.PhyloMaker")# add missing taxa to phylogeny backbone
library(ape) # phylogeny 
library(plyr) # data wrangling functions 
library(reshape2)# reshaping data into matrix
library(lubridate)# working with dates
library("FishLife") # traits data 
library("archetypes") # trait archetypes


# copy data from NARS IM
#########
# #NRSA0809-1819_siteinfo.tab
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809-1819_siteinfo.tab",
#           paste0(getwd(),"/Data/NRSA0809-1819_siteinfo.tab"))
# #NRSA0809-1819_fishCts_alltheNRSA.tab
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/NRSA0809-1819_fishCts_alltheNRSA.tab",
#           paste0(getwd(),"/Data/NRSA0809-1819_fishCts_alltheNRSA.tab"))
# #alltheNRSA_fish_taxa.tab
# file.copy("O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles/alltheNRSA_fish_taxa.tab",
#           paste0(getwd(), "/Data/alltheNRSA_fish_taxa.tab"))
################################################

SiteInfo <- read.table("Data/NRSA0809-1819_siteinfo.tab", sep="\t", header =T)

################################################################################
# Prepare taxonomy data file for all fish species collected from NRSA 
# Fish taxa unambiguously (e.g.NO hybrids or sp. or  were Amphibian) 
# identified to species, appended accepted name from FishBase using U.Taxonstand: 
# https://www.sciencedirect.com/science/article/pii/S2468265922000944
# https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13647

# analysis does not include Pomatomidae or PETROMYZONTIDAE
################################################################################

taxonomy <- read.table("Data/alltheNRSA_fish_taxa.tab", sep ="\t", header = T)
######
# roll up all unrecognized, but similar species to form species complexes
taxonomy[grep("SP. CF.", taxonomy$SPECIES, invert = F), "SPECIES"] <- 
  trimws(substring(taxonomy[grep("SP. CF.", taxonomy$SPECIES, invert = F),"SPECIES"],9))
taxonomy <- taxonomy[grep(" X ", taxonomy$GENUS, invert = T),]
taxonomy <- taxonomy[grep(" X |SP. | OR", taxonomy$SPECIES, invert = T),] 
taxonomy <- taxonomy[taxonomy$SPECIES != "SP.",]
taxonomy <- taxonomy[taxonomy$SPECIES != "SPECIES",]
taxonomy <- taxonomy[taxonomy$SPECIES != "HYBRID",]
taxonomy <- taxonomy[taxonomy$SPECIES != "", ] 
taxonomy <- taxonomy[taxonomy[, "HERP"] != "A", ]
taxonomy <- taxonomy[taxonomy$FAMILY != "Pomatomidae",] #Marine Pelagic Fish
taxonomy <- taxonomy[toupper(taxonomy$FAMILY) != "PETROMYZONTIDAE",] #Lamprey are not included in phylogeny
taxonomy$NRSA_SppName <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES)
taxonomy$NRSA_SppName <- toupper(taxonomy$NRSA_SppName)

# Roll up all subspecies to species level
x <- strsplit(taxonomy$NRSA_SppName, " ")
subspp <- taxonomy$NRSA_SppName[sapply(x, length) > 2]
taxonomy$NRSA_SppName <- unlist(lapply(x, function(i) paste(i[1], i[2])))
#######################################


################################################################################
# add accepted species names to autecology file 
# using U.Taxonstand r package. NRSA Species names may have 
# undergone a name change or spelling error. Any mismatched names 
# were checked manually
################################################################################

########
# database downloaded from package GitHub page 
FishBase <- read.xlsx("Data_Raw/Fishes_FishBase_database.xlsx")
matchName_Fishbase <- nameMatch(spList = taxonomy$NRSA_SppName, spSource = FishBase)
#matchName_Fishbase[matchName_Fishbase$Fuzzy,]

# Rename family to accepted family
names(matchName_Fishbase)[19] <- "Accepted_Family"

# Merged back to taxonomy file 
taxonomy <- merge(taxonomy, 
                  unique(matchName_Fishbase[,c("Submitted_Name", "Accepted_SPNAME", "Accepted_Family", "NOTE")]), 
                  by.x = "NRSA_SppName", by.y = "Submitted_Name", all.x = T)
taxonomy$Accepted_Genus <- sapply(strsplit(taxonomy$Accepted_SPNAME, " "), function(x) x[1])

# FNAME is the final name assigned for this study
taxonomy$FName <- taxonomy$Accepted_SPNAME
#########################################

# check instances where autecology names is different from Accepted_SPNAME 
# Confirm that accepted name is actually synonymous w/ autecology name. 
# goggled each Accepted_SPNAME, compared the common name. 
##############

updatedNames <- taxonomy[taxonomy$NRSA_SppName != toupper(taxonomy$Accepted_SPNAME),
                         c("TAXA_ID", "FINAL_NAME", "NRSA_SppName", 
                           "Accepted_Family", "Accepted_SPNAME", 
                           "FName", "NOTE")]

# For mismatched names, searched each record in AFS Names also results of google searches of 
# accepted names and nrsa names
updatedNames[order(updatedNames$FName),]

# The name given by NRSA was correct, as it matches AFS names
taxonomy$FName[taxonomy$NRSA_SppName == "ETHEOSTOMA MEADIAE"] <- "Etheostoma meadiae"
taxonomy$FName[taxonomy$NRSA_SppName == "MOXOSTOMA DUQUESNEI"] <- "Moxostoma Duquesnei"
taxonomy$FName[taxonomy$NRSA_SppName == "PANTOSTEUS LAHONTAN"] <- "Pantosteus lahontan"

# https://fieldguide.mt.gov/speciesDetail.aspx?elcode=AFC4E02440
taxonomy$FName[taxonomy$NRSA_SppName == "COTTUS BONDI"] <- "Cottus bairdii"
# https://explorer.natureserve.org/Taxon/ELEMENT_GLOBAL.2.102592/Etheostoma_spectabile
taxonomy$FName[taxonomy$NRSA_SppName == "ETHEOSTOMA PULCHELLUM"] <- "Etheostoma spectabile"

# older names match phylogeny/traits. Change FNAME to increase matching records (1:1 synonym)
# use old names for grouping: https://explorer.natureserve.org/Taxon/ELEMENT_GLOBAL.2.769114/Lepidomeda_copei
taxonomy$FName[taxonomy$FName == "Lepidomeda copei"] <- "Snyderichthys copei"
taxonomy$FName[taxonomy$FName == "Nothonotus aquali"] <- "Etheostoma aquali"
taxonomy$FName[taxonomy$FName == "Nothonotus sanguifluus"] <- "Etheostoma Sanguifluum"
############################################


################################################################################
# Search for instances in the taxonomy where separate NRSA taxa were combined 
# into a single taxa. The implication with combining records is that we will 
# be modeling a "complex" such that subspecies or individual (e.g. Cottus Bondi 
# and Cottus bairdii) taxa will not necessarily be distinguishable. All cutthroat
# trout are grouped Oncorhynchus clarkii, even though several subspecies have recently
# been elevated to species level
################################################################################

a <- split(taxonomy, taxonomy$FName)
a[unlist(lapply(a, function(s) nrow(s)>1))]

################################################################################
# Append occurrence data to FNames, 
################################################################################

# Compile site and species data downloaded from NARS IM
# hybrid and ambiguous taxa are removed from the analysis
# some taxa are listed in the taxonomy file but do not have 
# an occurrence record.
#########
# Master table
RawDataPath <- "O:/PRIV/CPHEA/PESD/COR/CORFILES/IM-TH007/data/im/allTheNRSA/data/tabfiles"
Fish_cnt <- RawDataPath |>
  paste0("/NRSA0809-1819_fishCts_alltheNRSA.tab")|>
  read.table(sep = "\t", header = T)

 
# fish count data by merging data with harmonized taxonomy file
Fish_cnt_Full <- merge(taxonomy[,c("TAXA_ID", "FName")], Fish_cnt, by = "TAXA_ID") |>
  mutate(YEAR = substring(DATE_COL, nchar(DATE_COL)-3, nchar(DATE_COL))) |>
  dplyr::select(c("UID","SITE_ID", "VISIT_NO", "YEAR", 
                  "UNIQUE_ID", "TAXA_ID", "FINAL_NAME", "FName", 
                  "TOTAL", "NON_NATIVE"))
#############################

#write.csv(Fish_cnt_Full, "Data/NRSAFish_Counts_AllYears_Complete.csv", row.names = F)


# check taxa that were collected by NRSA but not included in the analysis. These are taxa 
# that were either hybrids or ambiguous (e.g. UNKNOWN GENERA) 
removedTaxa <- names(sort(table(Fish_cnt[!Fish_cnt$TAXA_ID %in% taxonomy$TAXA_ID, "TAXA_ID"])))
taxonomy_check <- read.table("Data_Raw/allNRSA_FishTaxa.tab", sep ="\t", header = T)
taxonomy_check[taxonomy_check$TAXA_ID %in% removedTaxa,c("FINAL_NAME", "GENUS", "SPECIES")]


################################################################################
# Identify records for NRSA fish in megaphylogeny of fish
################################################################################

# using U.phylo pacaadd missing taxa to phylogenic tree
# https://www.sciencedirect.com/science/article/pii/S2468265922001329 to add
# taxa that were not included in the megephylogeny
# https://www.nature.com/articles/s41586-018-0273-1
# https://timetree.org/api/widget/citation_data/4453
##########
sp.list <- unique(Fish_cnt_Full$FName)
megatree <- read.tree("https://raw.githubusercontent.com/megatrees/fish_20221117/main/fish_megatree.tre")
gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/fish_20221117/main/fish_genus_list.csv', sep=",")
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
#######################################
result$sp.list$output.note
saveRDS(result, "Data/NRSAFish_Tree_Complete.rds")


# archetypes from FishLife as an alternative to FishMorph
# a more robust trait dataset that uses phylogenic imputation
#########
beta_iv = FishBase_and_Morphometrics$beta_gv[FishBase_and_Morphometrics$g_i,]
a <- archetypes(beta_iv, 3)
OPE <- setNames(data.frame(a$alphas, species = rownames(beta_iv)), c("O","P","E", "species"))
OPE$genus <- unlist(lapply(strsplit(OPE$species, " "), "[[",1))
OPE <- OPE[!duplicated(OPE$species),]

# create trait table
SpeciesList <- unique(Fish_cnt_Full$FName)
traits <- OPE[OPE$species%in%SpeciesList,]


# taxa that did not have a matching record
 add_taxa <- Fish_cnt_Full %>%
  select(FName) %>%
  distinct() %>%
  filter(!FName%in%OPE$species) %>%
  mutate(genus = unlist(lapply(strsplit(FName," "),"[[",1)))
  
  genMean <- sapply(unique(add_taxa$genus), 
                  simplify = F, function(x) 
                    colMeans(OPE[OPE$genus==x,
                                 c("O", "P", "E")]))
  genMean <- do.call(rbind, genMean)

  add_taxa <- merge(add_taxa, genMean, by.x = "genus", by.y = 0)

  # there was a name change in 
  add_taxa[add_taxa$FName == "Pantosteus lahontan", c("O", "P", "E")] <-
    OPE[OPE$species == "Catostomus platyrhynchus", c("O", "P", "E")]
  add_taxa[add_taxa$FName == "Snyderichthys copei", c("O", "P", "E")] <-
    OPE[OPE$FName == "Lepidomeda copei", c("O", "P", "E")]
  names(add_taxa)[2]<-"species"
##############################

traits <- rbind(traits, add_taxa)

write.csv(traits, "Data/NRSAFish_Traits_Complete.csv")


################################################################################
# Fish density data: Provided by Mike Mahon and Samantha Rumschlag (12/11/2023)
# update column names in fish density dataset to match FNames in taxonomy/count 
# dataset. 
################################################################################

FishCPUE <- read.csv("Data_Raw/NRSA_Fish_CPUE_MM_SR_12112023.csv")
###################
# merge tables by site c("SiteNumber","CollectionDate")
anyDuplicated(FishCPUE[, c("SiteNumber","CollectionDate")])

FishCPUEw_UID <- merge(SiteInfo[,c("UID", "DATE_COL", "UNIQUE_ID")],
                       FishCPUE,
                       by.x = c("UNIQUE_ID", "DATE_COL"),
                       by.y = c("SiteNumber", "CollectionDate"))

# all true -- no missing records
nrow(FishCPUEw_UID)[1] == nrow(FishCPUE)[1]

# species names in density dataset with >0 occurrence
# the dataset that was sent included species that were only collected from 
# USGS stations
dens.spp <- names(FishCPUEw_UID[,-c(1:22)])[colSums(FishCPUEw_UID[,-c(1:22)])>0]
FishCPUEw_UID <- FishCPUEw_UID[, c(names(FishCPUEw_UID[,c(1:22)]), dens.spp)]

# Update Phoxinus to Chrosomus and "Herichthys.cyanoguttatum" to 
# "Herichthys.cyanoguttatus"in density dataset to match NRSA FName, guidance 
# provided by MM 
# https://explorer.natureserve.org/Taxon/ELEMENT_GLOBAL.2.103219/Chrosomus_erythrogaster
names(FishCPUEw_UID)[grep(c("Phoxinus.erythrogaster|Phoxinus.oreas|Phoxinus.eos|Phoxinus.cumberlandensis|Phoxinus.neogaeus|Phoxinus.tennesseensis"), 
                          names(FishCPUEw_UID))]<- paste0("Chrosomus.", c("erythrogaster", "oreas", "eos", "cumberlandensis", "neogaeus", "tennesseensis"))
names(FishCPUEw_UID)[names(FishCPUEw_UID) == "Herichthys.cyanoguttatum"] <- "Herichthys.cyanoguttatus"

# In some instances, density names matched an unaccepted NRSA name 
# (i.e. NRSA_SppName). Updated column names in density dataset to match the FNAME,
# in taxonomy dataset

#melt to make it easier to update/group names
FishCPUE <- setNames(reshape2::melt(FishCPUEw_UID[,-c(1,2,4:22)], "UID"), 
                     c("UID","Name","Density"))
FishCPUE$Name <- gsub("\\."," ", FishCPUE$Name)

# Add FName column to FishCPU
FishCPUE$FName[FishCPUE$Name%in%unique(taxonomy$FName)] <- 
  FishCPUE$Name[FishCPUE$Name%in%unique(taxonomy$FName)]

# unique list of density taxa with unmatched record in NRSA taxonomy
fishwdenspp <- unique(as.character(FishCPUE$Name)[is.na(FishCPUE$FName)])

# Compare unmatched taxa to NRSA_SppName and update names in density dataset,
# Such that density data sp.names matches the NRSA FName
sppMatch <- taxonomy[taxonomy$NRSA_SppName %in% toupper(fishwdenspp), c("NRSA_SppName", "FName")]

for (i in 1:nrow(sppMatch)){
  # i <- 12
  FishCPUE$FName[toupper(FishCPUE$Name) %in% sppMatch$NRSA_SppName[i]] <- sppMatch$FName[i]    
}

#these missing taxa are associated with lamprey/rare taxa
unique(FishCPUE[is.na(FishCPUE$FName), "Name"])

FishCPUE <- FishCPUE[!is.na(FishCPUE$FName),]
#crosswalk between density data and occurrence data 
FNameCrosswalk <- unique(FishCPUE[,c("Name","FName")])


taxonomy$Density[taxonomy$FName %in% unique(FishCPUE$FName)] <- "Y"
taxonomy$Density[!taxonomy$FName %in% unique(FishCPUE$FName)] <- "N"
#################################################
write.csv(FishCPUE, "Data_Raw/CleanData/NRSAFish_Density_AllYears.csv", row.names = F)

#crosswalk sent to MM and SR to compare with density dataset
write.csv(FNameCrosswalk, "Data_Raw/CleanData/FNameCrosswalk.csv", row.names=F)
write.csv(taxonomy, "Data_Raw/CleanData/NRSAFish_Taxonomy.csv", row.names = F)
