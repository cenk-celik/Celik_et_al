home <- system("echo ${HOME}/", intern = TRUE)
main <- paste0(home,"nmrc/")

#----Load libraries & tweak processors----
library(Seurat)

seurat <- readRDS('seurat.rds')

#----Split by samples----
split.list <- SplitObject(seurat, split.by = "samples")

barcodes <- list()

for (i in 1:length(split.list)) {
  condition <- names(split.list[i])
  seurat <- unlist(split.list[[i]])
  barcodes[[i]] <- colnames(split.list[[i]])
  
  write.table(unlist(barcodes[[i]]),
            file = paste0(names(split.list[i]), '_barcodes.csv'),
            col.names = F,
            row.names = F
            )
  
}

sessionInfo()
