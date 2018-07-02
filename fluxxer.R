#!/usr/bin/env Rscript

#load necessary libraries
suppressMessages(library(rsalvador))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input CSV file name", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file prefix", metavar="character")
  )
  
usage_string = paste(
  "fluxxer.R -i input.csv -o output_prefix\n\n",
  "Input CSV must be in a tidy format with one plate count per line and the column ", 
  "headings: strain, plate, fraction, and CFU\n\n",
  "The columns must contain:\n",
  "  strain: name of the strain tested\n",
  "  plate: the type of plate count, either {selective|s} or {nonselective|count|ns}\n",
  "  fraction: the fraction of the culture that was plated\n",
  "  CFU: the number of colonies counted on the plate\n\n",
  "Fraction should generally be 1 for the selective plates (if you plated the entire culture).\n\n",
  "Fraction for the nonselective (count) plates should be equal to the ratio of the volume plated (P) ",
  "to the culture volume (C) divided by the dilution factor (D), or Fraction = P/(C*D).\n\n",
  "As an example, if you had 200 µl cultures and put the entire volume into 10 ml of saline in a ",
  "first dilution tube, then transferred 1 µl from this tube to another 10 ml in a second dilution ",
  "tube, and finally plated and counted the number of cells in 50 µl of this second dilution tube: ",
  "P = 50 µl, C = 200 µL, and D = (200 µl / 10,000 µl) * (1 µl / 10,000 µl) = 2E-6. ",
  "Therefore, Fraction = 50 µl / (200 µL * 2E-6) = 5E–7",
  sep = ""
) 

opt_parser = OptionParser(usage=usage_string, option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Supply -i argument for input CSV file", call.=FALSE)
}

calculateMutRate <- function(filename, output_prefix)
{
  
  #read in file specified. Must be in same directory
  #for testing
  #data <- read_csv("example_dataset_2.csv") 
  data = read_csv(filename)  
  
  #do some checks of the input files to expand abbreviations
  data$plate = tolower(data$plate)
  data = mutate(data, plate = ifelse( (plate == "n") | (plate == "ns") | (plate == "count"), "nonselective", plate))
  data = mutate(data, plate = ifelse(plate == "s", "selective", plate))
  
  data$strain = as.factor(data$strain)
  data$plate = as.factor(data$plate)

  strains = levels(data$strain)
  
  #identify # of strains, use to build empty data frame
  num_strains <- length(strains)
  cat("Found", num_strains, "strains:\n")
  cat(paste(strains, sep=", "), "\n\n")
  
  output_data <- tibble()
  
  #cycle through each column to calculate mutatation rate and confidence  
  for(this.strain in strains) {
    cat("\nSTRAIN:", this.strain, "\n")
    #locate Non_selective separator
    this.strain.data = data %>% filter(strain==this.strain)

    #extract selective values
    selective.rows = this.strain.data %>% filter(plate=="selective")
    nonselective.rows = this.strain.data %>% filter(plate=="nonselective")
    num_selective = nrow(selective.rows)
    num_nonselective = nrow(nonselective.rows)
    
    cat("Number of selective plate counts:", num_selective, "\n") 
    cat("Number of nonselective plate counts:", num_nonselective, "\n")
    
    
    if (num_selective == 0 || num_nonselective == 0 ) {
      cat("***ERROR! Did not find plate counts for selective/nonselective. Skipping strain.\n")
    }
    
    nonselective_cell_counts = mean(nonselective.rows$CFU/nonselective.rows$fraction)
    cat("Estimated cells per culture:", nonselective_cell_counts, "(", nrow(nonselective.rows), "nonselective plates )\n")
    
    #all selective plates must have the same fraction
    selective_fraction_list = selective.rows %>% count(fraction)
    if (nrow(selective_fraction_list) > 1) {
      cat("***ERROR! Multiple fractions found for selective plates. Skipping strain.\n")
      next
    }
    selective_fraction = selective_fraction_list$fraction[1]
    cat("Fraction or efficiency of selective cultures plated (e):", selective_fraction, "\n")
    
    if (selective_fraction == 1) {
      m = newton.LD(selective.rows$CFU)
    } else {
      m = newton.LD.plating(selective.rows$CFU, e=selective_fraction)
    }
    
    mu = m / nonselective_cell_counts
    cat("Maximum likelihood mutation rate (mu):", mu, "\n")
    
    if (selective_fraction == 1) {
      CI = confint.LD(selective.rows$CFU, alpha=0.05)/nonselective_cell_counts
    } else {
      CI = confint.LD.plating(selective.rows$CFU, alpha=0.05, e=selective_fraction)/nonselective_cell_counts
    }
    cat("         95% confidence interval (mu): [", CI[1], ",", CI[2] , "]\n")
    
    output_data = rbind(output_data, data.frame(strain = this.strain, num_nonselective = num_nonselective, num_selective = num_selective, selective_fraction = selective_fraction, mu = mu, CI.95.lower = CI[1], CI.95.higher = CI[2]))

  }
  
  write_csv(output_data, paste0(output_prefix,".output.csv"))
  ##make chart for pretty values
  plot <- ggplot(output_data, aes(x = strain, y = mu)) +
    geom_point() +
    geom_linerange(aes(ymin = CI.95.lower, ymax = CI.95.higher)) +
    scale_y_log10() + 
    ggtitle("Mutation Rates") + 
    xlab("Strains") +
    ylab("Mutation rate MLE") +
    annotation_logticks(sides = "l")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  save_plot(paste0(output_prefix, ".plot.pdf"), plot)
  
}

calculateMutRate(filename = opt$input, output_prefix = opt$output)
