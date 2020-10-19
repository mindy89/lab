#read file
#remember to change the pathway of file = "..." depend on different device

#note that it is conserved sites of miRNA on Target Scan which was written by order, not download data
temp <- read.csv(file = "C:/Users/ªô¨K©û/Documents/Lab/data/survival/KRAS/PAAD/temp.txt", header = F, sep = "\n", stringsAsFactor = F)
#clinical data download from TCGA
clinical <- read.csv(file = "C:/Users/ªô¨K©û/Documents/Lab/data/survival/KRAS/PAAD/clinical.tsv", header = T, sep = "\t", stringsAsFactor = F)
#98 miRNA from KRAS gene
for(l in 1:98){
  miR <- temp[l, 1]
  #samplesheet download from TCGA
  samplesheet <- read.csv(file = "C:/Users/ªô¨K©û/Documents/Lab/data/survival/KRAS/PAAD/gdc_sample_sheet.2020-09-11.tsv", header = T, sep = "\t", stringsAsFactor = F)
  
  for(i in 1:260){ #cases from clinical raw data
    for(j in 1:135){ #cases from samplesheet
      if(i %% 2 == 1){ #clinical raw data repeat twice, thus divided by 2
        if(clinical[i, 2] == samplesheet[j, 6]){ #if find same caseID
          #get each cases' miRNA expression
          mir <- read.csv(file = paste0("C:/Users/ªô¨K©û/Documents/Lab/data/survival/KRAS/PAAD/", samplesheet[j, 1], "/", samplesheet[j, 2]), header = T, sep = "\t")
          #total of 1881 miRNA
          for(k in 1:1881){
            if(mir[k, 1] == miR){
              find.row <- k
              break
            }
          }
          samplesheet[j, 5] <- mir[find.row, 2] #miRNA expression
          samplesheet[j, 4] <- clinical[i, 16] #vital status
          if(samplesheet[j, 4] == "Alive"){
            samplesheet[j, 4] <- 0
            samplesheet[j, 3] <- clinical[i, 48] #days to last follow up
          }
          else if(samplesheet[j, 4] == "Dead"){
            samplesheet[j, 4] <- 1
            samplesheet[j, 3] <- clinical[i, 10] #days to death
          }
        }
      }
    }
  }
  
  #if vital status = "0"|"--", delete that row
  #note that if the situation occurred continuously, create another loop or repeat running the function
  #134 stands for the cases left
  for(i in 1:134){
    if(samplesheet[i, 3] == "0"){
      samplesheet <- samplesheet[-i, ]
      }
  else if(samplesheet[i, 3] == "'--"){
      samplesheet <- samplesheet[-i, ]
    }
  }
  
  
  #output
  
  #rename the name of column
  colnames(samplesheet) <- c("ID", "file.name", "days_to_last_follow_up", "vital_status", miR)
  #adjusting the column
  samplesheet[, 1] <- samplesheet[, 6]
  #removing useless columns
  for(i in 1:3){
    samplesheet[, 6] <- NULL
  }
  #output .csv
  write.csv(samplesheet, file = paste0("C:/Users/ªô¨K©û/Documents/Lab/data/KRAS/PAAD/", miR, ".csv"))
  
  #read file
  PAAD <- read.csv(file = paste0("C:/Users/ªô¨K©û/Documents/Lab/data/survival/KRAS/PAAD/", miR, ".csv"), stringsAsFactor = F)
  #removing useless column
  PAAD <- PAAD[, -1]
  #get the median of miRNA expression
  M <- median(PAAD[, 5])
  #separate miRNA expression into high and low expression group
  for(i in 1:134){
    if(PAAD[i, 5] > M){
      PAAD[i, 6] <- ("high")
    }
    else{
      PAAD[i, 6] <- ("low") 
    }
  }
  #rename the name of column
  colnames(PAAD) <- c("ID", "file.name", "days_to_last_follow_up", "vital_status", miR, "expression")
  #install.packages("ggplot2")
  #install.packages("survminer")
  #include lib
  library("survminer")
  require("survival")
  fit <- survfit(Surv(PAAD$days_to_last_follow_up, PAAD$vital_status) ~ PAAD$expression)
  #to draw survival curve and view in plots
  ggsurvplot(fit, data = PAAD, xlim = c(0, 2300), xlab = "days to last follow up", pval = TRUE)
  output <- ggsurvplot(fit, data = PAAD, xlim = c(0, 2300), xlab = "days to last follow up", pval = TRUE, print = TRUE)
  ggsave(filename = paste0("C:/Users/ªô¨K©û/Documents/Lab/data/KRAS/PAAD/PAAD_", miR, ".png"), plot = print(output), width = 5, height = 4)
}

#user interaction. Use while running specific miRNA, not all. Remember not to run miR <- temp[l, 1]
miR <- readline(prompt = "miRNA:")