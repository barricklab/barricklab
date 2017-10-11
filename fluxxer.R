#load necessary libraries
library(rsalvador)
library(tidyverse)
library(cowplot)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file prefix", metavar="character"),
  make_option(c("-p", "--volumeplated"), type = "integer", default = 50, 
                help="volume (in microliters) plated, default is 50"),
  make_option(c("-v", "--originalvolume"), type = "integer", default = 200, 
                help="original volume (in microliters), default is 200"),
  make_option(c("-d", "--dilutionfactor"), type = "double", default = 2e-6, 
                help="dilution factor, default is 2e-6")

  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Supply -i argument for input file", call.=FALSE)
}

#confidence interval function
confind <- function(x, y){
  S <- as.vector(x)
  NS <- as.vector(y)
  return(confint.LD(S, alpha=0.05)/NS[1])
}

#golden.benchmark.LD function for calculating mutation rate
golden <- function(x, y){
  #S <- as.vector(x)
  #NS <- as.vector(y)
  return(golden.benchmark.LD(x, y, tol = 1e-09, max.iter = 100, show.iter = FALSE))
}


calculateMutRate <- function(filename, vol_plated = 50, dilution_factor = 2e-6, original_volume = 200, out_pfx ){


  #read in file specified. Must be in same directory
  #for testing
  #data <- read_csv("example_dataset_2.csv") 
  data <- read_csv(filename)
  #extract sample names
  strains <- colnames(data)
  
  #identify # of strains, use to build empty data frame
  num_strains <- length(strains)
  data_output <- tibble(strain = rep("", num_strains), MLE = rep(0, num_strains), confint.lower = rep(0, num_strains), confint.higher = rep(0, num_strains))

#cycle through each column to calculate mutatation rate and confidence  
for(i in 1:length(strains)){
  #locate Non_selective separator
  sep <- which(data[,i] == "Non_Selective")
  
  #extract selective values
  selective <- as.numeric(data[1:sep-1, i][[1]])
  selective <- selective[!is.na(selective)]
  #extract non selective values
  raw_non_selective <- as.numeric(data[-(1:sep), i][[1]])
  per_non_selective <- ((mean(raw_non_selective, na.rm = TRUE) / vol_plated) / dilution_factor ) * original_volume
  
  #calculate non selective based on dilution factor (placeholder till correct calculation implemented)
  #match length of nonselective to selective values
  non_selective <- rep(per_non_selective, length(selective))
  print(selective)
  print(non_selective)
  print(strains[i])
  values <- c(strains[i], golden(selective, non_selective), confind(selective, non_selective))
  data_output[i,] <- values
}
##write data frame to output csv
data_output$MLE <- as.numeric(data_output$MLE)
data_output$confint.lower <- as.numeric(data_output$confint.lower)
data_output$confint.higher <- as.numeric(data_output$confint.higher)
data_output$strain <- factor(data_output$strain, strains)


print(data_output)
write_csv(data_output, paste0(out_pfx,"_output.csv"))
##make chart for pretty values
plot <- ggplot(data_output, aes(x = strain, y = MLE)) +
  geom_point() +
  geom_linerange(aes(ymin = confint.lower, ymax = confint.higher)) +
  scale_y_log10() + 
  ggtitle("Mutation Rates") + 
  xlab("Strains") +
  ylab("Mutation rate MLE") +
  annotation_logticks(sides = "l")

save_plot(paste0(out_pfx, "_chart.pdf"), plot)

}

calculateMutRate(filename = opt$input, vol_plated = opt$volumeplated, dilution_factor = opt$dilutionfactor, original_volume = opt$originalvolume, out_pfx = opt$out)
