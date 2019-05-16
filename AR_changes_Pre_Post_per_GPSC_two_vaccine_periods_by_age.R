#Calculate the changes in prevlance of antibiotic resistance between two vaccine periods with age group stratification
library(epiR)
library(pwr)

# Read and subset dataset
Metadata <- read.delim2("GPS_Brazil_metadata.csv", header = T, sep = ',') #input the metadata in csv format
AR_props <- subset(Metadata, select=c(GPSC, Vaccine_Period__autocolour, Age_group__autocolour, WGS_PEN_SIR_Meningitis__autocolour, WGS_CHL_SIR__autocolour, WGS_COT_SIR__autocolour, WGS_ERY_SIR__autocolour, WGS_TET_SIR__autocolour))

                                            
# Create an output dataframe
AR_Pre_PostPCV <- data.frame("Age" =0, "Antibiotic"=0, "GPSC"=0, "Pre-PCV_R" =0, "Pre-PCV_total" =0, "Post-PCV_R" =0, "Post-PCV_total"=0, "p"=0)
Pre_R <- 0
Pre_S <- 0
Post_R <-0
Post_S <-0
Pre_R_proportion <- 0
Pre_S_proportion <- 0
Post_R_proportion <- 0
Post_S_proportion <- 0
non_susceptible <- c("I","R")

#Power test and Fisher's exact test
GPSCs <- unique(AR_props$GPSC) 
age <- unique (AR_props$Age_group__autocolour)

for (k in age){
  AR_props_age <- subset(AR_props, Age_group__autocolour == k)

  for(i in GPSCs){
    AR_props_age_GPSC <- subset(AR_props_age, GPSC==i )
  
    for( j in grep("WGS_", colnames(AR_props_age_GPSC), value = T))
    {
    
        Pre_Post_AR <- AR_props_age_GPSC[c(j, "Vaccine_Period__autocolour")]
        pcv <- unique(AR_props_age_GPSC$Vaccine_Period__autocolour)
        
        Pre_R <- nrow(subset(Pre_Post_AR,Vaccine_Period__autocolour=="Pre-PCV" & Pre_Post_AR[,1] %in% c("I", "R")))
        Pre_S <- nrow(subset(Pre_Post_AR,Vaccine_Period__autocolour=="Pre-PCV" & Pre_Post_AR[,1] =="S"))
        Post_R <- nrow(subset(Pre_Post_AR,Vaccine_Period__autocolour=="Post-PCV10" & Pre_Post_AR[,1] %in% c("I", "R")))
        Post_S <- nrow(subset(Pre_Post_AR,Vaccine_Period__autocolour=="Post-PCV10" & Pre_Post_AR[,1] =="S"))
        
        Pre_R_proportion <- Pre_R / (Pre_R + Pre_S + 0.001) + 0.001
        Pre_S_proportion <- Pre_S / (Pre_R + Pre_S + 0.001) + 0.001
        Post_R_proportion <- Post_R / (Post_R + Post_S + 0.001) + 0.001
        Post_S_proportion <- Post_S / (Post_R + Post_S + 0.001) + 0.001
      
        AR <- matrix(c(Pre_R,Pre_S,Post_R, Post_S), nrow = 2, byrow = T)
        AR_proportion <- matrix(c(Pre_R_proportion,Pre_S_proportion,Post_R_proportion,Post_S_proportion), nrow = 2, byrow = T)
        pwr <- pwr.chisq.test(w = ES.w2(AR_proportion), power = 0.8, df = 1, sig.level = 0.05)
        
        if(sum(AR) > pwr$N) 
        {
        rownames(AR) <- pcv
        colsum_AR <- colSums(AR)
        row <- AR[1,, drop=F]
        by2 <- matrix(c(row[1,1],row[1,2],colsum_AR[1]-row[1,1],colsum_AR[2]-row[1,2]),ncol=2,byrow=TRUE)
        p <- fisher.test(by2, alternative)$p.value
        AR_Pre_PostPCV <- rbind(AR_Pre_PostPCV, c(k,j, i, Pre_R,Pre_R+Pre_S,Post_R,Post_R+Post_S,p))
        }
    
    }
  }
}

#Stratification by age and write output in csv
AR_Pre_PostPCV_less_than_5 <- subset(AR_Pre_PostPCV, Age=="<5")
Z <- nrow(AR_Pre_PostPCV_less_than_5)
AR_Pre_PostPCV_less_than_5$newcolumn <- p.adjust(AR_Pre_PostPCV_less_than_5$p, "BH", n= Z)
names(AR_Pre_PostPCV_less_than_5)[names(AR_Pre_PostPCV_less_than_5) == "newcolumn"] <- "adjusted_p"
write.csv(AR_Pre_PostPCV_less_than_5, file="AR_Pre_PostPCV_less_than_5.csv", row.names = FALSE)

AR_Pre_PostPCV_more_than_5 <- subset(AR_Pre_PostPCV, Age==">=5")
Y <- nrow(AR_Pre_PostPCV_more_than_5)
AR_Pre_PostPCV_more_than_5$newcolumn <- p.adjust(AR_Pre_PostPCV_more_than_5$p, "BH", n= Y)
names(AR_Pre_PostPCV_more_than_5)[names(AR_Pre_PostPCV_more_than_5) == "newcolumn"] <- "adjusted_p"
write.csv(AR_Pre_PostPCV_more_than_5, file="AR_Pre_PostPCV_more_than_5.csv", row.names = FALSE)



