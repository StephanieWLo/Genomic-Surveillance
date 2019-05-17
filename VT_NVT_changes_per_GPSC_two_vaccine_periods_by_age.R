#Detect proportional changes in GPSC between vaccine periods among vaccine serotype (VT) and non-vaccine serotype (NVT) separately
#Tab3_NVT_prop_perGPSC.R from Rebecca Gladstone https://github.com/rgladstone detect such proportional changes too. This script added a Power test to calculate the minimum sample size for comparisons and age stratification.

library(epiR)
library(plyr)
library(pwr)

#Prepare input dataframe
Metadata <- read.delim2("GPS_Brazil_metadata.csv", header = T, sep = ',')
data <- subset(Metadata, select=c("GPSC", grep("Vaccine_Period", colnames(Metadata), value = TRUE), grep("PCV13_Status", colnames(Metadata), value = TRUE), grep("^Age_group", colnames(Metadata), value = TRUE)))

# Create an output dataframe
GPSC_Pre_PostPCV <- data.frame("Age"=0, "Vaccine_Status" =0, "GPSC"=0, "Pre-PCV_GPSC" =0, "Pre-PCV_total" =0, "Post-PCV_GPSC" =0, "Post-PCV_total"=0, "p"=0)
Pre_GPSC <- 0
Pre_other <- 0
Post_GPSC <-0
Post_other <-0
Pre_GPSC_proportion <- 0
Pre_other_proportion <- 0
Post_GPSC_proportion <- 0
Post_other_proportion <- 0


#Power test and Fisher's exact test
age <- unique(data$Age_group__autocolour)
vaccine_status <- unique(data$PCV13_Status__autocolour)

for (k in age){
  data_age <- subset(data, data$Age_group__autocolour == k)
    for (j in vaccine_status){
      data_age_status <- subset(data_age, data_age$PCV13_Status__autocolour == j)
      GPSCs <- unique(data_age_status$GPSC)
        for (i in GPSCs){
          data_age_status_GPSC <- subset(data_age_status, GPSC == i)
          pcv <- sort(unique(data_age_status_GPSC$Vaccine_Period__autocolour))
          Pre_GPSC <- nrow(subset(data_age_status_GPSC, Vaccine_Period__autocolour==pcv[1][1]))
          Post_GPSC <- nrow(subset(data_age_status_GPSC, Vaccine_Period__autocolour==pcv[2][1]))
          Pre_other <- nrow(subset(data_age_status,GPSC !=i & Vaccine_Period__autocolour == pcv[1][1]))
          Post_other <- nrow(subset(data_age_status,GPSC !=i & Vaccine_Period__autocolour == pcv[2][1]))
          
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
            GPSC_Pre_PostPCV <- rbind(GPSC_Pre_PostPCV, c(k, j, i, Pre_GPSC,Pre_GPSC+Pre_other,Post_GPSC,Post_GPSC+Post_other,p))
            
          }
        }  
    }
}


#Stratification by age and write output in csv
GPSC_Pre_PostPCV_less_than_5 <- subset(GPSC_Pre_PostPCV, Age == "<5")
GPSC_Pre_PostPCV_less_than_5$newcolumn <- p.adjust(GPSC_Pre_PostPCV_less_than_5$p, "BH", n = nrow(GPSC_Pre_PostPCV_less_than_5))
names(GPSC_Pre_PostPCV_less_than_5)[names(GPSC_Pre_PostPCV_less_than_5) == "newcolumn"] <- "adjusted_p"
write.csv(GPSC_Pre_PostPCV_less_than_5, file ="GPSC_Pre_PostPCV_less_than_5.csv", row.names = FALSE)

GPSC_Pre_PostPCV_more_than_5 <- subset(GPSC_Pre_PostPCV, Age == ">=5")
GPSC_Pre_PostPCV_more_than_5$newcolumn <- p.adjust(GPSC_Pre_PostPCV_more_than_5$p, "BH", n = nrow(GPSC_Pre_PostPCV_more_than_5))
names(GPSC_Pre_PostPCV_more_than_5)[names(GPSC_Pre_PostPCV_more_than_5) == "newcolumn"] <- "adjusted_p"
write.csv(GPSC_Pre_PostPCV_more_than_5, file ="GPSC_Pre_PostPCV_more_than_5.csv", row.names = FALSE)
