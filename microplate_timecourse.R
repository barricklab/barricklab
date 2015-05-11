library(ggplot2)

##
## REQUIRED INPUT FILES
##
## input/OD600.tab
## input/GFP.tab
## input/RFP.tab
##
## Each table for a different signal (OD600, GFP, RFP) must be put in
## a different tab-delimited file for input. The example_input directory 
## contains an example of creating these files from an Excel spreadsheet.
##
## input/key.tab 
## 
## This file must be created by the user. It include information
## about the sample in each well. Currently the number of replicate
## background wells for an experiment must be equal to or greater than
## the number of replicates for all conditions, and the number of wells
## for each condition must be equal.

input_path = "input"
output_path = "output"

plate_key_file_path = file.path(input_path, "key.tab")
plate_OD600_file_path = file.path(input_path, "OD600.tab")
plate_RFP_file_path = file.path(input_path, "RFP.tab")
plate_GFP_file_path = file.path(input_path, "GFP.tab")

##ggplot2 defaults
theme_set(theme_bw(base_size = 24))
theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
line_thickness = 0.8

#################################

key.table = read.table(plate_key_file_path, header=T)
OD600.table = read.table(plate_OD600_file_path, header=F, row.names=1)
RFP.table = read.table(plate_RFP_file_path, header=F, row.names=1)
GFP.table = read.table(plate_GFP_file_path, header=F, row.names=1)

#this makes a list of all levels that are not none in aaRS
aaRS.list = levels(droplevels(subset(key.table, aaRS != "none")$aaRS))


for(this.aaRS in aaRS.list) {
  
  mega.table = data.frame()
  
  this.aa = as.character(subset(key.table, (aaRS==this.aaRS) & (AA != "none"))$AA[1])
  
  #Find the wells that are just amino acid and LB
  background.wells = subset(key.table, (aaRS=="none") & (AA == this.aa) )
  
  ## take the first three
  background.noAA.wells = subset(key.table, (medium=="LB") & (AA=="none") & (aaRS=="none") )
 # background.noAA.wells = background.noAA.wells[1:3,]
  
  #How many clones?
  clone.list = levels(droplevels(subset(key.table, (aaRS==this.aaRS) & (clone != "none"))$clone))

  
  
  for(this.clone in clone.list) {
  
    pFRYC.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA!="none") & (plasmid=="pFRYC"))
    pFRY.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA!="none") & (plasmid=="pFRY"))
    
    # note: this code will die if there are not exactly the right
    # number of wells for pFRYC, pFRY, and the background
    this.bg.OD600.table = OD600.table[as.character(background.wells$well),]
    this.bg.GFP.table = GFP.table[as.character(background.wells$well),]
    this.bg.RFP.table = RFP.table[as.character(background.wells$well),]
    
    this.pFRYC.OD600.table = OD600.table[as.character(pFRYC.wells$well),]
    this.pFRYC.GFP.table = GFP.table[as.character(pFRYC.wells$well),]
    this.pFRYC.RFP.table = RFP.table[as.character(pFRYC.wells$well),]
    
    this.pFRY.OD600.table = OD600.table[as.character(pFRY.wells$well),]
    this.pFRY.RFP.table = RFP.table[as.character(pFRY.wells$well),]
    this.pFRY.GFP.table = GFP.table[as.character(pFRY.wells$well),]    
    
    pFRYC.noAA.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA=="none") & (plasmid=="pFRYC"))
    pFRY.noAA.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA=="none") & (plasmid=="pFRY"))
    
    this.bg.OD600.noAA.table = OD600.table[as.character(background.noAA.wells$well),]
    this.bg.GFP.noAA.table = GFP.table[as.character(background.noAA.wells$well),]
    this.bg.RFP.noAA.table = RFP.table[as.character(background.noAA.wells$well),]
    
    this.pFRYC.OD600.noAA.table = OD600.table[as.character(pFRYC.noAA.wells$well),]
    this.pFRYC.GFP.noAA.table = GFP.table[as.character(pFRYC.noAA.wells$well),]
    this.pFRYC.RFP.noAA.table = RFP.table[as.character(pFRYC.noAA.wells$well),]
    
    this.pFRY.OD600.noAA.table = OD600.table[as.character(pFRY.noAA.wells$well),]
    this.pFRY.RFP.noAA.table = RFP.table[as.character(pFRY.noAA.wells$well),]
    this.pFRY.GFP.noAA.table = GFP.table[as.character(pFRY.noAA.wells$well),]   
    
    #create new mega.table rows for all signals minus background
    for (i in 1:nrow(this.pFRY.OD600.table)) {
    
      for (j in 1:ncol(OD600.table)) {
      
        new.row = data.frame(
          aaRS = this.aaRS,
          clone = this.clone,
          replicate = i,
          time=OD600.table["Time",j],
          OD600.bg = this.bg.OD600.table[i,j],
          GFP.bg = this.bg.GFP.table[i,j],
          RFP.bg = this.bg.RFP.table[i,j],
          OD600.pFRYC = this.pFRYC.OD600.table[i,j],
          RFP.pFRYC = this.pFRYC.RFP.table[i,j],
          GFP.pFRYC = this.pFRYC.GFP.table[i,j],
          OD600.pFRY = this.pFRY.OD600.table[i,j],
          RFP.pFRY = this.pFRY.RFP.table[i,j],
          GFP.pFRY = this.pFRY.GFP.table[i,j],
          
          OD600.bg.noAA = this.bg.OD600.noAA.table[i,j],
          GFP.bg.noAA = this.bg.GFP.noAA.table[i,j],
          RFP.bg.noAA = this.bg.RFP.noAA.table[i,j],
          OD600.pFRYC.noAA = this.pFRYC.OD600.noAA.table[i,j],
          RFP.pFRYC.noAA = this.pFRYC.RFP.noAA.table[i,j],
          GFP.pFRYC.noAA = this.pFRYC.GFP.noAA.table[i,j],
          OD600.pFRY.noAA = this.pFRY.OD600.noAA.table[i,j],
          RFP.pFRY.noAA = this.pFRY.RFP.noAA.table[i,j],
          GFP.pFRY.noAA = this.pFRY.GFP.noAA.table[i,j]        
          
        )
      
        mega.table = rbind(mega.table, new.row)
      }
    }
  }
  
  #subtract backgrounds
  
  mega.table$OD600.pFRYC.net = mega.table$OD600.pFRYC - mega.table$OD600.bg
  mega.table$RFP.pFRYC.net = mega.table$RFP.pFRYC - mega.table$RFP.bg
  mega.table$GFP.pFRYC.net = mega.table$GFP.pFRYC - mega.table$GFP.bg
  mega.table$OD600.pFRY.net = mega.table$OD600.pFRY - mega.table$OD600.bg
  mega.table$RFP.pFRY.net = mega.table$RFP.pFRY - mega.table$RFP.bg
  mega.table$GFP.pFRY.net = mega.table$GFP.pFRY - mega.table$GFP.bg
  
  mega.table$RFP.pFRYC.normalized = mega.table$RFP.pFRYC.net / mega.table$OD600.pFRYC
  mega.table$GFP.pFRYC.normalized = mega.table$GFP.pFRYC.net / mega.table$OD600.pFRYC
  mega.table$pFRYC.ratio = mega.table$RFP.pFRYC.normalized / mega.table$GFP.pFRYC.normalized
  
  mega.table$RFP.pFRY.normalized = mega.table$RFP.pFRY.net / mega.table$OD600.pFRY
  mega.table$GFP.pFRY.normalized = mega.table$GFP.pFRY.net / mega.table$OD600.pFRY
  mega.table$pFRY.ratio = mega.table$RFP.pFRY.normalized / mega.table$GFP.pFRY.normalized
  
  mega.table$decoding.efficiency = mega.table$pFRY.ratio / mega.table$pFRYC.ratio
  
  mega.table$OD600.pFRYC.noAA.net = mega.table$OD600.pFRYC.noAA - mega.table$OD600.bg.noAA
  mega.table$RFP.pFRYC.noAA.net = mega.table$RFP.pFRYC.noAA - mega.table$RFP.bg.noAA
  mega.table$GFP.pFRYC.noAA.net = mega.table$GFP.pFRYC.noAA - mega.table$GFP.bg.noAA
  mega.table$OD600.pFRY.noAA.net = mega.table$OD600.pFRY.noAA - mega.table$OD600.bg.noAA
  mega.table$RFP.pFRY.noAA.net = mega.table$RFP.pFRY.noAA - mega.table$RFP.bg.noAA
  mega.table$GFP.pFRY.noAA.net = mega.table$GFP.pFRY.noAA - mega.table$GFP.bg.noAA

  mega.table$RFP.pFRYC.noAA.normalized = mega.table$RFP.pFRYC.noAA.net / mega.table$OD600.pFRYC.noAA
  mega.table$GFP.pFRYC.noAA.normalized = mega.table$GFP.pFRYC.noAA.net / mega.table$OD600.pFRYC.noAA
  mega.table$pFRYC.noAA.ratio = mega.table$GFP.pFRYC.noAA.normalized / mega.table$RFP.pFRYC.noAA.normalized
  
  mega.table$RFP.pFRY.noAA.normalized = mega.table$RFP.pFRY.noAA.net / mega.table$OD600.pFRY.noAA
  mega.table$GFP.pFRY.noAA.normalized = mega.table$GFP.pFRY.noAA.net / mega.table$OD600.pFRY.noAA
  mega.table$pFRY.noAA.ratio = mega.table$GFP.pFRY.noAA.normalized / mega.table$RFP.pFRY.noAA.normalized
  
  mega.table$decoding.efficiency.noAA = mega.table$pFRY.noAA.ratio / mega.table$pFRYC.noAA.ratio
  
  mega.table$misincorporation.value =  mega.table$decoding.efficiency.noAA  /  mega.table$decoding.efficiency
  
  ### ADD CODE HERE
  
  #mainly for debugging
  write.csv(mega.table, file.path(output_path, paste(this.aaRS, "_mega_table.csv", sep="")))
         
  mega.table$replicate = factor(mega.table$replicate)
  mega.table$clone_replicate = paste(mega.table$clone, mega.table$replicate, sep="_")
  
  p = ggplot(mega.table, aes(x=time, y=mega.table$decoding.efficiency, color=clone, linetype=replicate))
  p + geom_line(shape=0, size=line_thickness) + coord_cartesian(ylim=c(-1,3))
  ggsave(filename=file.path(output_path, paste(this.aaRS, "_decoding.efficiency.AA.pdf", sep="")))
 
 p = ggplot(mega.table, aes(x=time, y=mega.table$decoding.efficiency.noAA, color=clone, linetype=replicate))
 p + geom_line(shape=0, size=line_thickness) + coord_cartesian(ylim=c(-1,3))
 ggsave(filename=file.path(output_path, paste(this.aaRS,"_decoding.efficiency.noAA.pdf", sep="")))
 
 p = ggplot(mega.table, aes(x=time, y=mega.table$misincorporation.value, color=clone, linetype=replicate))
 p + geom_line(shape=0, size=line_thickness) + coord_cartesian(ylim=c(-1,3))
 ggsave(filename=file.path(output_path, paste(this.aaRS,"_misincorporation.value.noAA.pdf", sep="")))
}



