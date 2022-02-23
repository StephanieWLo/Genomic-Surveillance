#GPS_QC
#pf assembly -t file -i lanes.txt -s -> lanes.txt.assemblyfind_stats.csv
#pf qc -t file -i lanes.txt -L S -s -T -> qc_summary.csv
#pf data -t file -i lanes.txt -s

#Get the following csv files using pf commands
assembly <- read.delim2("lanes.txt.assemblyfind_stats.csv", header = T, sep = ',')
kraken <- read.delim2("qc_summary.csv", header = T, sep = ',')
raw_reads <- read.delim2("lanes.txt.pathfind_stats.csv", header = T, sep = ',')
hetsnp <- read.delim2("combinedNewStats.txt", header = F, sep = "\t")

colnames(hetsnp) <- c("lane_id", "hetsnp_count", "homo_count", "ratio")

temp1 = merge(assembly, kraken, by.x = "Lane", by.y = "Species", all.x = TRUE)
temp2 = merge(temp1, raw_reads, by.x = "Lane", by.y = "Lane.Name", all.x = TRUE)
temp3 = merge(temp2, hetsnp, by.x ="Lane", by.y ="lane_id", all.x = TRUE)

data <- subset(temp3, select=c("Lane", "Streptococcus.pneumoniae", "Total.Length", "No.Contigs", "Genome.Covered", "Depth.of.Coverage", "hetsnp_count"))

for (i in unique(colnames(data))[-1]){
  data[,i] <- as.numeric(as.character(data[,i]))
}

data$qc <- with(data, ifelse(data$Streptococcus.pneumoniae >60 &
                               data$Total.Length < 2300000 &
                               data$Total.Length > 1900000 &
                               data$No.Contigs < 500 &
                               data$Genome.Covered >60 &
                               data$Depth.of.Coverage >20 &
                               data$hetsnp_count <220,
                                "Pass", "Fail"))

pass <- subset(data, data$qc == "Pass")
pass$qc <- with(pass, ifelse(pass$No.Contigs <100, "Pass Plus", "Pass"))
fail <- subset(data, data$qc == "Fail")
qc_output <- rbind(pass, fail)
write.csv(qc_output, file = "QC_output.csv", row.names = FALSE)
