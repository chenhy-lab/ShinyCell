library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(gridExtra) 
library(RColorBrewer)
library(readr)
# library(ShinyCell)
library(reticulate)
use_virtualenv("/ocean/home/chenhy/project/Shiny/ShinyCell/src/my_shiny_py")
# python virtualenv.pyz -p /ocean/apps/miniconda3/envs/py37R4/bin/python my_shiny_py

# use_condaenv("/ocean/apps/miniconda3/envs/py37R4")
# use_python(python='/ocean/apps/miniconda3/envs/py37R4/bin/python')
py_config()
py_module_available("anndata")


# h5ad/loom
seu <- '/ocean/home/wangyaqiong/project/single_cell/mouse/embryonic/processed/Eye_adata_scale.h5ad'

data_path <- "/data/chenhy/Shiny/ShinyCell/"
# ProjectName <- "test"
ProjectName <- "mouse_test"
file_list <- list.files('/ocean/home/chenhy/project/Shiny/ShinyCell/src/R', pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list)source(file_path)
scConf <- createConfig(seu)
makeShinyApp(seu, scConf, gene.mapping = T,shiny.dir = sprintf("%s/%s",data_path,ProjectName))






# presto-test
# RDS-Seurat
setwd('~/project/Shiny/ShinyCell/')
seu <- readRDS('/ocean/home/chenhy/project/dongxr/daidan/data/sce.combineExpr_T_cell.rds')
file_list <- list.files('/ocean/home/chenhy/project/Shiny/ShinyCell/src/R/V1.1', pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {source(file_path)}
scConf <- createConfig(seu)
makeShinyApp(seu, scConf, gene.mapping = T,shiny.dir = 'data/project/test/test/test')

setwd('~/project/Shiny/ShinyCell/data/project/test/test/test')
library(presto)
library(rhdf5)

library(dplyr)



sampling_proportions <- c(1, 0.5, 0.05, 0.005)
names(sampling_proportions) <- c('1-1000', '1001-10000', '10001-100000', '100001-1000000')
select_sampling_proportions <- function(cell_counts) {
  proportions <- rep(0, length(cell_counts))
  for (i in 1:length(cell_counts)) {
    count <- cell_counts[i]
    for (range in sampling_ranges) {
      range_bounds <- strsplit(range, '-')[[1]]
      lower_bound <- as.numeric(range_bounds[1])
      upper_bound <- as.numeric(range_bounds[2])
      if (count >= lower_bound && count <= upper_bound) {
        proportions[i] <- sampling_proportions[range]
        break
      }
    }
  }
  return(proportions)
}

cell_counts <- sc1meta %>%
  dplyr::group_by(.data[[group]]) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::pull(count)

# Selection of sampling ratios based on cell number intervals
use_sampling_proportions <- select_sampling_proportions(cell_counts)

# Grouping at each level and selecting appropriate sampling ratios for sampling based on cell count intervals
set.seed(156468)
sampled_metadata <- sc1meta %>%
  dplyr::group_by(.data[[group]]) %>%
  dplyr::sample_frac(use_sampling_proportions, replace = FALSE) %>% 
  dplyr::ungroup()
cell_name_index <- match(sampled_metadata$sampleID, sc1meta$sampleID)


X_matrix <- Matrix::Matrix(h5data$read(args = list(quote(expr=),cell_name_index))) 



X_matrix <- Seurat::GetAssayData(seu, assay = "RNA", slot = "data")
group_by <- "orig.ident"
if (is.null(group_by)) {
  y <- Seurat::Idents(seu)
} else {
  y <- Seurat::FetchData(seu, group_by) %>% unlist %>% as.character()
}



sc1meta <- readRDS('data/project/test/test/test/sc1meta.rds')
sc1gene = readRDS("data/project/test/test/test/sc1gene.rds")




res <- FindAllMarkers(metadata = sc1meta,group_by = 'RNA_snn_res.0.2')


res <- FindMarkers(metadata = sc1meta,group_by = 'orig.ident',group.1='case',group.2='control')

res2 <- TopMarker(res)
