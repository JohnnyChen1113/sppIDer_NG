#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
require(ggplot2)

args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("No strain name provided as an argument")
}
strainName <- args[1]

# Define paths and read data
workingDir <- file.path(getwd(), "")
file_base <- file.path(workingDir, strainName)
MQ_file <- paste0(file_base, "_MQ.txt")
chiSqOutName <- paste0(file_base, "_MQ_chiSquared.txt")
summaryOutputName <- paste0(file_base, "_MQsummary.txt")
plotOutputName <- paste0(file_base, "_plotMQ.pdf")

# Helper functions
read_and_prepare_data <- function(file_path) {
  MQdf <- read.table(file_path, header=TRUE)
  MQdf$Species <- factor(MQdf$Species, levels = unique(MQdf$Species))
  MQdf$Species[MQdf$Species == "*"] <- "Unmapped"
  return(MQdf)
}

calculate_statistics <- function(data) {
  maps <- subset(data, Species != "Unmapped")
  MQscores <- subset(data, MQscore > 0)
  total_reads <- sum(data$count)
  mapped_reads <- sum(maps$count)
  prop_zero <- (sum(subset(data, MQscore == 0)$count) / total_reads) * 100
  avg_MQ <- weighted.mean(data$MQscore, data$count)
  median_MQ <- median(rep(data$MQscore, times=data$count))
  
  return(list(total_reads = total_reads, mapped_reads = mapped_reads, prop_zero = prop_zero,
              avg_MQ = avg_MQ, median_MQ = median_MQ))
}

write_summary <- function(summary_data, file_name) {
  write(paste(strainName, "Num reads = ", toString(summary_data$total_reads), sep=""), file=file_name)
  write(paste(strainName, "Num mapped reads = ", toString(summary_data$mapped_reads), sep=""), file=file_name, append=TRUE)
  write(paste(strainName, "Unmapped reads = ", toString(round(summary_data$prop_zero, digits=2)), "%", sep=""), file=file_name, append=TRUE)
  write(paste(strainName, "Average MQ = ", toString(summary_data$avg_MQ), sep=""), file=file_name, append=TRUE)
  write(paste(strainName, "Median MQ = ", toString(summary_data$median_MQ), sep=""), file=file_name, append=TRUE)
}

# Main data processing
MQdf <- read_and_prepare_data(MQ_file)
stats <- calculate_statistics(MQdf)
write_summary(stats, summaryOutputName)

# More complex analysis and plots can be added here
# ...

# Set up plot parameters and create plots
pdf(plotOutputName)
ggplot(MQdf, aes(x=factor(Species), fill=Species)) +
  geom_bar(aes(y=..prop..), stat="count") +
  labs(title=paste(strainName, "Mapping Quality Distribution"), x="Species", y="Proportion") +
  theme_minimal()
dev.off()
