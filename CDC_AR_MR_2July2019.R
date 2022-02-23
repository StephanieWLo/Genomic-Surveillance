#Format CDC SPN AR output into Microreact metadata format
library(data.table)
library (tidyverse)

# Read and subset dataset
AR <- read.delim2("TABLE_SPN_new_Typing_Results.csv", header = TRUE, sep = ",") #input the metadata in csv format
AR <- AR %>% rename(lane_id = Sample)
data <- subset(AR, select=c("lane_id", "PBP1A","PBP2B","PBP2X",grep("WGS_", colnames(AR), value = TRUE), "EC", "Cot", "Tet", "FQ","Other"))

#Concatenate pbp profiles and autocolour them
data$PBP1A_2B_2X__autocolour <- paste(data$PBP1A, data$PBP2B, data$PBP2X, sep ="--")

#Colour columns
colnames(data)[which(names(data) == "Pili")] <- "Pili__autocolour"
#colnames(data)[which(names(data) == "Cot")] <- "Cot__autocolour"
colnames(data)[which(names(data) == "Tet")] <- "Tet__autocolour"
colnames(data)[which(names(data) == "FQ")] <- "FQ__autocolour"
#colnames(data)[which(names(data) == "Other")] <- "Other__autocolour"

#Concatenate the SIGN and SIR
for (k in grep("_SIGN", colnames(data), value = TRUE)) 
{
  data[, gsub("_SIGN", "", k)] <- paste(gsub("^=$","", data[,k]), data[,gsub("_SIGN", "", k)])
}
  
#delete SIGN columns
drops <- grep("_SIGN", colnames(data), value = TRUE)
data <- data[ , !names(data) %in% drops]

drops <- grep("WGS_Serotype", colnames(data), value = TRUE)
data <- data[ , !names(data) %in% drops]

#Add Microreact colour
for (i in grep("_SIR", colnames(data), value = TRUE)) 
{
  data[,paste0(i, "__colour")] <- "transparent"
  data[which(data[,i] %in% "R"), paste0(i, "__colour")] <- "#ff2722"
  data[which(data[,i] %in% "S"), paste0(i, "__colour")] <- "#0069ec"
  data[which(data[,i] %in% "I"), paste0(i, "__colour")] <- "#f797b1"
}

#Extract resistance genotypes
data[,paste0("ermB")] <- "neg"
ermB_row <- which(data$EC %like% "ERM")
data$ermB[ermB_row] <- "pos"

data[,paste0("ermB__colour")] <- "transparent"
ermB_colour_row <- which(data$ermB == "pos")
data$ermB__colour[ermB_colour_row] <- "#ff2722"
ermB_colour_row <- which(data$ermB == "neg")
data$ermB__colour[ermB_colour_row] <- "#0069ec"

data[,paste0("mefA")] <- "neg"
mefA_row <- which(data$EC %like% "MEF")
data$mefA[mefA_row] <- "pos"

data[,paste0("mefA__colour")] <- "transparent"
mefA_colour_row <- which(data$mefA == "pos")
data$mefA__colour[mefA_colour_row] <- "#ff2722"
mefA_colour_row <- which(data$mefA == "neg")
data$mefA__colour[mefA_colour_row] <- "#0069ec"

data[,paste0("folA_I100L")] <- "neg"
folA_I100L_row <- which(data$Cot %like% "I20L")
data$folA_I100L[folA_I100L_row] <- "pos"

data[,paste0("folA_I100L__colour")] <- "transparent"
folA_I100L_colour_row <- which(data$folA_I100L == "pos")
data$folA_I100L__colour[folA_I100L_colour_row] <- "#ff2722"
folA_I100L_colour_row <- which(data$folA_I100L == "neg")
data$folA_I100L__colour[folA_I100L_colour_row] <- "#0069ec"

data[,paste0("folP")] <- "neg"
folP_row <- grep("FOLP", data$Cot)
data$folP[folP_row] <- as.character(data$Cot[folP_row])
data$folP <- gsub(".*FOLP", "FOLP", data$folP)
colnames(data)[which(names(data) == "folP")] <- "folP__autocolour"

data[,paste0("cat")] <- "neg"
cat_row <- grep("CAT", data$Other)
data$cat[cat_row] <- "pos"

data[,paste0("cat__colour")] <- "transparent"
cat_colour_row <- which(data$cat == "pos")
data$cat__colour[cat_colour_row] <- "#ff2722"
cat_colour_row <- which(data$cat == "neg")
data$cat__colour[cat_colour_row] <- "#0069ec"

drops <- c("WGS_AMP", "WGS_AMP_SIR", "WGS_AMP_SIR__colour", "WGS_CPT", "WGS_CPT_SIR", "WGS_CPT_SIR__colour", "WGS_ZOX", "WGS_ZOX_SIR", "WGS_ZOX_SIR__colour", "WGS_FOX", "WGS_FOX_SIR", "WGS_FOX_SIR__colour", "WGS_CIP", "WGS_CIP_SIR", "WGS_CIP_SIR__colour", "WGS_DAP", "WGS_DAP_SIR", "WGS_DAP_SIR__colour")
trim_data <- data[ ,!(names(data) %in% drops)]

#data <- data[,c(1:5,62,6:61,63:88)]

write.csv(trim_data, file ="UK_AMR.csv", row.names = FALSE)

#Merge AR output with Metadata file
#Metadata <- read.delim2("Bangladesh_metadata_9Aug2019_1a.csv", header = T, sep = ',')
#AR <- read.delim2("Bangladesh_AR_output_formatted.csv", header = T, sep = ',')
#Final = merge(Metadata, AR, by.x = "Sample", by.y = "Sample", all.x = TRUE)
#write.csv(Final, file = "Russia_metadata_8Sep2019.csv", row.names = FALSE)

