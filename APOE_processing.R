#this creates a column with both APOE alleles combined
apoeres <- read.csv("~/APOERES.csv") 
apoeres$apoe <- paste(paste("E", apoeres$APGEN1, sep = ""), paste("E", apoeres$APGEN2, sep = ""), sep = "/")
