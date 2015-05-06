library(ggplot2)

plate_key_file_path = "key.tab"
plate_OD600_file_path = "OD600.tab"
plate_RFS_file_path = "RFS.tab"
plate_GFS_file_path = "GFS.tab" 

##ggplot2 defaults
theme_set(theme_bw(base_size = 24))
theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
line_thickness = 0.8

#################################

key.table = read.table(plate_key_file_path, header=T)
OD600.table = read.table(plate_OD600_file_path, header=F, row.names=1 )
RFS.table = read.table(plate_RFS_file_path, header=F, row.names=1)
GFS.table = read.table(plate_GFS_file_path, header=F, row.names=1)

#this makes a list of all levels that are not none in aaRS
aaRS.list = levels(droplevels(subset(key.table, aaRS != "none")$aaRS))


for(this.aaRS in aaRS.list) {
  
  mega.table = data.frame()
  
  this.aa = as.character(subset(key.table, (aaRS==this.aaRS) & (AA != "none"))$AA[1])
  
  #Find the wells that are just amino acid and LB
  background.wells = subset(key.table, (aaRS=="none") & (AA == this.aa) )
  
  ## take the first three
  background.noAA.wells = subset(key.table, (medium=="LB") & (AA=="none") & (aaRS=="none") )
  background.noAA.wells = background.noAA.wells[1:3,]
  
  #How many clones?
  clone.list = levels(droplevels(subset(key.table, (aaRS==this.aaRS) & (clone != "none"))$clone))

  
  
  for(this.clone in clone.list) {
  
    pFRYC.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA!="none") & (plasmid=="pFRYC"))
    pFRY.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA!="none") & (plasmid=="pFRY"))
    
    # note: this code will die if there are not exactly the right
    # number of wells for pFRYC, pFRY, and the background
    this.bg.OD600.table = OD600.table[as.character(background.wells$well),]
    this.bg.GFS.table = GFS.table[as.character(background.wells$well),]
    this.bg.RFS.table = RFS.table[as.character(background.wells$well),]
    
    this.pFRYC.OD600.table = OD600.table[as.character(pFRYC.wells$well),]
    this.pFRYC.GFS.table = GFS.table[as.character(pFRYC.wells$well),]
    this.pFRYC.RFS.table = RFS.table[as.character(pFRYC.wells$well),]
    
    this.pFRY.OD600.table = OD600.table[as.character(pFRY.wells$well),]
    this.pFRY.RFS.table = RFS.table[as.character(pFRY.wells$well),]
    this.pFRY.GFS.table = GFS.table[as.character(pFRY.wells$well),]    
    
    pFRYC.noAA.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA=="none") & (plasmid=="pFRYC"))
    pFRY.noAA.wells = subset(key.table, (aaRS==this.aaRS) & (clone==this.clone) & (AA=="none") & (plasmid=="pFRY"))
    
    this.bg.OD600.noAA.table = OD600.table[as.character(background.noAA.wells$well),]
    this.bg.GFS.noAA.table = GFS.table[as.character(background.noAA.wells$well),]
    this.bg.RFS.noAA.table = RFS.table[as.character(background.noAA.wells$well),]
    
    this.pFRYC.OD600.noAA.table = OD600.table[as.character(pFRYC.noAA.wells$well),]
    this.pFRYC.GFS.noAA.table = GFS.table[as.character(pFRYC.noAA.wells$well),]
    this.pFRYC.RFS.noAA.table = RFS.table[as.character(pFRYC.noAA.wells$well),]
    
    this.pFRY.OD600.noAA.table = OD600.table[as.character(pFRY.noAA.wells$well),]
    this.pFRY.RFS.noAA.table = RFS.table[as.character(pFRY.noAA.wells$well),]
    this.pFRY.GFS.noAA.table = GFS.table[as.character(pFRY.noAA.wells$well),]   
    
    #create new mega.table rows for all signals minus background
    for (i in 1:nrow(this.bg.OD600.table)) {
    
      for (j in 1:ncol(OD600.table)) {
      
        new.row = data.frame(
          aaRS = this.aaRS,
          clone = this.clone,
          replicate = i,
          time=OD600.table["Time",j],
          OD600.bg = this.bg.OD600.table[i,j],
          GFS.bg = this.bg.GFS.table[i,j],
          RFS.bg = this.bg.RFS.table[i,j],
          OD600.pFRYC = this.pFRYC.OD600.table[i,j],
          RFS.pFRYC = this.pFRYC.RFS.table[i,j],
          GFS.pFRYC = this.pFRYC.GFS.table[i,j],
          OD600.pFRY = this.pFRY.OD600.table[i,j],
          RFS.pFRY = this.pFRY.RFS.table[i,j],
          GFS.pFRY = this.pFRY.GFS.table[i,j],
          
          OD600.bg.noAA = this.bg.OD600.noAA.table[i,j],
          GFS.bg.noAA = this.bg.GFS.noAA.table[i,j],
          RFS.bg.noAA = this.bg.RFS.noAA.table[i,j],
          OD600.pFRYC.noAA = this.pFRYC.OD600.noAA.table[i,j],
          RFS.pFRYC.noAA = this.pFRYC.RFS.noAA.table[i,j],
          GFS.pFRYC.noAA = this.pFRYC.GFS.noAA.table[i,j],
          OD600.pFRY.noAA = this.pFRY.OD600.noAA.table[i,j],
          RFS.pFRY.noAA = this.pFRY.RFS.noAA.table[i,j],
          GFS.pFRY.noAA = this.pFRY.GFS.noAA.table[i,j]        
          
        )
      
        mega.table = rbind(mega.table, new.row)
      }
    }
  }
  
  #subtract backgrounds
  
  mega.table$OD600.pFRYC.net = mega.table$OD600.pFRYC - mega.table$OD600.bg
  mega.table$RFS.pFRYC.net = mega.table$RFS.pFRYC - mega.table$RFS.bg
  mega.table$GFS.pFRYC.net = mega.table$GFS.pFRYC - mega.table$GFS.bg
  mega.table$OD600.pFRY.net = mega.table$OD600.pFRY - mega.table$OD600.bg
  mega.table$RFS.pFRY.net = mega.table$RFS.pFRY - mega.table$RFS.bg
  mega.table$GFS.pFRY.net = mega.table$GFS.pFRY - mega.table$GFS.bg
  
  mega.table$RFS.pFRYC.normalized = mega.table$RFS.pFRYC.net / mega.table$OD600.pFRYC
  mega.table$GFS.pFRYC.normalized = mega.table$GFS.pFRYC.net / mega.table$OD600.pFRYC
  mega.table$pFRYC.ratio = mega.table$RFS.pFRYC.normalized / mega.table$GFS.pFRYC.normalized
  
  mega.table$RFS.pFRY.normalized = mega.table$RFS.pFRY.net / mega.table$OD600.pFRY
  mega.table$GFS.pFRY.normalized = mega.table$GFS.pFRY.net / mega.table$OD600.pFRY
  mega.table$pFRY.ratio = mega.table$RFS.pFRY.normalized / mega.table$GFS.pFRY.normalized
  
  mega.table$incorporation.value = mega.table$pFRY.ratio / mega.table$pFRYC.ratio
  
  mega.table$OD600.pFRYC.noAA.net = mega.table$OD600.pFRYC.noAA - mega.table$OD600.bg.noAA
  mega.table$RFS.pFRYC.noAA.net = mega.table$RFS.pFRYC.noAA - mega.table$RFS.bg.noAA
  mega.table$GFS.pFRYC.noAA.net = mega.table$GFS.pFRYC.noAA - mega.table$GFS.bg.noAA
  mega.table$OD600.pFRY.noAA.net = mega.table$OD600.pFRY.noAA - mega.table$OD600.bg.noAA
  mega.table$RFS.pFRY.noAA.net = mega.table$RFS.pFRY.noAA - mega.table$RFS.bg.noAA
  mega.table$GFS.pFRY.noAA.net = mega.table$GFS.pFRY.noAA - mega.table$GFS.bg.noAA

  mega.table$RFS.pFRYC.noAA.normalized = mega.table$RFS.pFRYC.noAA.net / mega.table$OD600.pFRYC.noAA
  mega.table$GFS.pFRYC.noAA.normalized = mega.table$GFS.pFRYC.noAA.net / mega.table$OD600.pFRYC.noAA
  mega.table$pFRYC.noAA.ratio = mega.table$RFS.pFRYC.noAA.normalized / mega.table$GFS.pFRYC.noAA.normalized
  
  mega.table$RFS.pFRY.noAA.normalized = mega.table$RFS.pFRY.noAA.net / mega.table$OD600.pFRY.noAA
  mega.table$GFS.pFRY.noAA.normalized = mega.table$GFS.pFRY.noAA.net / mega.table$OD600.pFRY.noAA
  mega.table$pFRY.noAA.ratio = mega.table$RFS.pFRY.noAA.normalized / mega.table$GFS.pFRY.noAA.normalized
  
  mega.table$incorporation.value.noAA = mega.table$pFRY.noAA.ratio / mega.table$pFRYC.noAA.ratio
  
  mega.table$incorporation.value.ratio =  mega.table$incorporation.value /  mega.table$incorporation.value.noAA
  
  ### ADD CODE HERE
  
  #mainly for debugging
  write.csv(mega.table, paste(this.aaRS, "_mega_table.csv", sep=""))
         
  mega.table$replicate = factor(mega.table$replicate)
  mega.table$clone_replicate = paste(mega.table$clone, mega.table$replicate, sep="_")
  
  p = ggplot(mega.table, aes(x=time, y=mega.table$GFS.pFRYC.noAA.normalized, color=clone, linetype=replicate))
  p + geom_line(shape=0, size=line_thickness) + coord_cartesian()
  ggsave(filename=file.path(paste(this.aaRS, ".pdf", sep="")))
}



