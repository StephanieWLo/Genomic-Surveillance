#Detect proportional changes in GPSC between vaccine periods among vaccine serotype (VT) and non-vaccine serotype (NVT) separately
#Tab3_NVT_prop_perGPSC.R from Rebecca Gladstone https://github.com/rgladstone detect such proportional changes too. This script added a Power test to calculate the minimum sample size for comparisons.

library(epiR)
library(plyr)
library(pwr)

#Prepare input dataframe
Metadata <- read.delim2("GPS_Argentina_metadata.csv", header = T, sep = ',')
data <- subset(Metadata, select=c("GPSC", "Vaccine_Period", "Vaccine_Status")) 

# Create an output dataframe
GPSC_Pre_PostPCV <- data.frame("Vaccine_Period" =0, "Vaccine_Status" =0, "GPSC"=0, "Pre-PCV_GPSC" =0, "Pre-PCV_total" =0, "Post-PCV_GPSC" =0, "Post-PCV_total"=0, "p"=0)
#VT_Pre_PostPCV <- data.frame("Vaccine" =0, "GPSC"=0, "Pre-PCV_GPSC_VT" =0, "Pre-PCV_total_VT" =0, "Post-PCV_GPSC_VT" =0, "Post-PCV_total_VT"=0, "p")
Pre_GPSC <- 0
Pre_other <- 0
Post_GPSC <-0
Post_other <-0
Pre_GPSC_proportion <- 0
Pre_other_proportion <- 0
Post_GPSC_proportion <- 0
Post_other_proportion <- 0


#Power test and Fisher's exact test
GPSCs <- unique(data$GPSC) 
vaccine_period <- unique(data$Vaccine_Period)
vaccine_status <- unique(data$Vaccine_Status)

for (k in vaccine_status){
  data_status <- subset(data, data$Vaccine_Status == k)
    for (j in vaccine_period){
      data_status_period <- subset(data_status, data_status$Vaccine_Period !=j)
        for (i in GPSCs){
          data_status_period_GPSC <- subset(data_status_period, data_status_period$GPSC == i)
          pcv <- unique(data_status_period_GPSC$Vaccine_Period)
          Pre_GPSC <- nrow(subset(data_status_period_GPSC, Vaccine_Period==pcv[1][1]))
          Post_GPSC <- nrow(subset(data_status_period_GPSC, Vaccine_Period==pcv[2][1]))
          Pre_other <- nrow(subset(data_status_period,GPSC !=i & Vaccine_Period == pcv[1][1]))
          Post_other <- nrow(subset(data_status_period,GPSC !=i & Vaccine_Period == pcv[2][1]))
          
          Pre_GPSC_proportion <- Pre_GPSC / (Pre_GPSC + Pre_other + 0.001) + 0.001
          Post_GPSC_proportion <- Post_GPSC / (Post_GPSC + Post_other + 0.001) + 0.001
          Pre_other_proportion <- Pre_other / (Pre_GPSC + Pre_other + 0.001) + 0.001
          Post_other_proportion <- Post_other / (Post_GPSC + Post_other + 0.001) + 0.001
          
          GPSC <- matrix(c(Post_GPSC,Pre_GPSC,Post_other, Pre_other), nrow = 2, byrow = T)
          GPSC_proportion <- matrix(c(Post_GPSC_proportion,Pre_GPSC_proportion,Post_other_proportion, Pre_other_proportion), nrow = 2, byrow = T)
          pwr <- pwr.chisq.test(w = ES.w2(GPSC_proportion), power = 0.8, df = 1, sig.level = 0.05)
          
          if(sum(GPSC) > pwr$N) 
          {
            p <- fisher.test(GPSC, alternative)$p.value
            pcv_compare <- unique(data_status_period_GPSC$Vaccine_Period)
            GPSC_Pre_PostPCV <- rbind(GPSC_Pre_PostPCV, c(toString(pcv_compare),k, i, Pre_GPSC,Pre_GPSC+Pre_other,Post_GPSC,Post_GPSC+Post_other,p))
            
          }
        }  
    }
}
GPSC_Pre_PostPCV <- GPSC_Pre_PostPCV[-1,]

#Stratification by vaccine_status and write output in csv
GPSC_VT_Pre_PostPCV <- subset(GPSC_Pre_PostPCV, Vaccine_Status == "PCV")
GPSC_VT_Pre_PostPCV$newcolumn <- p.adjust(GPSC_VT_Pre_PostPCV$p, "BH", n = nrow(GPSC_VT_Pre_PostPCV))
names(GPSC_VT_Pre_PostPCV)[names(GPSC_VT_Pre_PostPCV) == "newcolumn"] <- "adjusted_p"
write.csv(GPSC_VT_Pre_PostPCV, file ="GPSC_VT_Pre_PostPCV.csv", row.names = FALSE)

GPSC_NVT_Pre_PostPCV <- subset(GPSC_Pre_PostPCV, Vaccine_Status == "NVT")
GPSC_NVT_Pre_PostPCV$newcolumn <- p.adjust(GPSC_NVT_Pre_PostPCV$p, "BH", n = nrow(GPSC_NVT_Pre_PostPCV))
names(GPSC_NVT_Pre_PostPCV)[names(GPSC_NVT_Pre_PostPCV) == "newcolumn"] <- "adjusted_p"
write.csv(GPSC_NVT_Pre_PostPCV, file ="GPSC_NVT_Pre_PostPCV.csv", row.names = FALSE)