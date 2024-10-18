BACON
===
## About BACON

BACteria COmmunicatioN (BACON), a method that utilizes single-microbe RNA sequencing data to infer, visualize and analyze quorum sensing on single species and human gut microbiome.<br>
![](https://i-blog.csdnimg.cn/direct/caba390ea99344999e624e0aa05cba35.jpeg#pic_center)



## Installation
You can install the latest version of BACON like so:

```r
devtools::install_github("qu-wx/BACON")
```

## Usage

```r
library(BACON)
## create BACON object
example_BACON_object <- createBACON(object = as.matrix(gene_count_matrix),group.by = metadata$select_vector,database = database)
## calculation of quorum sensing networks
example_BACON_object <- run_BACON(object = example_BACON_object,M=100)
## aggregating the network
qs_aggregated <- net_aggregation(example_BACON_object@net)
## visulization
netVisual_circle_BACON(example_BACON_object@net$select_qs)
```
## Tutorials
[Database construction for cross-species](https://github.com/qu-wx/BACON-database/tree/main) 
[Inference and visualization of quorum sensing network](https://github.com/qu-wx/BACON-database/blob/main/tutor/Inference%20and%20visualization%20of%20quorum%20sensing%20network.md)
## R dependencies (tested and recommended)

```r
R >= 4.1.0  
dplyr >= 1.0.9
data.table >= 1.14.2  
Seurat >= 4.1.0  
SeuratObject >= 4.1.0  
igraph >= 1.3.4  
ggplot2 >= 3.3.6  
ComplexHeatmap >= 2.8.0  
circlize >= 0.4.14      
ggalluvial >= 0.12.3  
```
![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/2255f5911f2c409581f98c7226873a46.jpeg#pic_center)


