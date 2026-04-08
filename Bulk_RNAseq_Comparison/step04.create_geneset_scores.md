---
title: "hCBO: geneset stats "
author: "Qiang (Danny) Chen"
date: "2025-09-19"
output:
  html_document:
    keep_md: true
    highlight: tango
    code_folding: hide
    number_sections: yes
    theme: united
    toc: yes
    toc_float:
      collapsed: yes
      smooth_scroll: no
  word_document:
    toc: yes
  pdf_document:
    toc: yes
  bibliography: biblio.bib
vignette: |
  %\VignetteIndexEntry{RNASeq} %\VignetteEngine{knitr::rmarkdown} %\VignetteEncoding{UTF-8} %\usepackage[utf8]{inputenc}
---

We integrated SCTL samples at day60 in our 2 month-old hCBO data with the following two datasets and identified DEGs among cerebellar markers across protocols. 

For the quadrato lab, we focused on three samples below (they generated three replicates at day 60):  
[GSE247974](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247974)
GSM7904010: organoid D1, scRNA-seq  
GSM7904011: organoid D2, scRNA-seq  
GSM7904012: organoid D3, scRNA-seq  

For the pasca lab, we only focused on one sample (at day 72) below:  
[GSE233574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233574)
GSM7430706: OrganoidScreen; sample11; FGF2-50  

There is only one sample (GSM7430706) in GSE233574. We did sub-sampling in cells of GSM7430706 and calculated psudo-bulk of sub-samples for comparisons. 



``` r
# knitr::opts_chunk$set(tidy=FALSE, cache=FALSE, echo=TRUE, dev="png", message=FALSE, error=FALSE, warning=FALSE)
knitr::opts_chunk$set(echo=TRUE, warning=F, message=F)

library("tidyr")
library("dplyr")
library("ggplot2")
library("knitr")
library("pheatmap")
library("DT")
library("car")
library("lme4")
library("emmeans")

rm(list=ls())

d_meta = "/data/qiangchen/Projects/NCATS/SCTL/SeungMiRyu/hCBO/Data/"
d_proj = "/data/qiangchen/Projects/NCATS/SCTL/SeungMiRyu/hCBO/SCTLday60_GSE247974_GSE233574/"
d_dat  = paste0(d_proj, "Data/")
d_res  = paste0(d_proj, "Results/")
```

# Load data


``` r
f_cb = paste0(d_meta, "Full_Brain_Region_Marker_Table_Ryu_et_al.csv")
gene_data = read.csv(f_cb)

f_sample = paste0(d_dat, "sample_info.SCTLday60_GSE247974_GSE233574.tsv")
sampleInfo = read.table(f_sample, sep="\t", header=T)

sctl_samples = sampleInfo[sampleInfo$Protocol=="SCTLday60", "ID"]

f_norm = paste0(d_dat, "norm_counts.SCTLday60_GSE247974_GSE233574.CerebellarGenes.rds")
norm_counts = as.data.frame(readRDS(f_norm))

norm_data = as.data.frame(cbind(Gene=rownames(norm_counts), norm_counts))

norm_data = pivot_longer(norm_data, 
    cols = 2:ncol(norm_data), 
    names_to = "ID", 
    values_to = "Count")

norm_data = merge(norm_data, sampleInfo, by="ID")
```

# Geneset stats for cell types:

Bergmann glia  
Deep cerebellar nuclei (excitatory/inhibitory)  
Endothelial  
Golgi cell  
Granule cell / lineage  
Microglia  
Molecular layer interneuron (basket/stellate)  
Oligodendrocyte / OPC  
Pericyte / SMC   
Purkinje neuron  
Unipolar brush cell (UBC)  


Mixed effect model was used to compare protocols


``` r

celltypes = c("Bergmann glia", "Deep cerebellar nuclei (excitatory/inhibitory)", 
    "Endothelial", "Golgi cell", "Granule cell / lineage", "Microglia", 
    "Molecular layer interneuron (basket/stellate)", "Oligodendrocyte / OPC", 
    "Pericyte / SMC", "Purkinje neuron", "Unipolar brush cell (UBC)")

test_gene_data = gene_data[gene_data$CellType_subclass %in% celltypes, ]

test_count_data = norm_data[norm_data$Gene %in% test_gene_data$MarkerGene, ]
dim(test_count_data)
## [1] 990   7

contrasts = list()
contrasts[[1]] = c("SCTLday60", "GSE233574")
contrasts[[2]] = c("SCTLday60", "GSE247974")

nconts = length(contrasts)

statList = lapply(1:nconts, function(icont) {
    
    print(icont)

    statList = lapply(celltypes, function(celltype) {

        print(celltype)
    
        geneset = unique(test_gene_data[test_gene_data$CellType_subclass==celltype, "MarkerGene"])

        cur_protocols = as.vector(contrasts[[icont]])
        
        cur_counts = test_count_data[which( (test_count_data$Gene %in% geneset) & (test_count_data$Protocol %in% cur_protocols) ), ]
        cur_counts$Protocol = factor(cur_counts$Protocol, levels=cur_protocols)
        cur_counts$Gene = as.factor(cur_counts$Gene)

        model = lmer(Count ~ Protocol + (1 | Gene), data = cur_counts)

        emm = emmeans(model, ~ Protocol) 
        sum_cont = summary(contrast(emm, list(sample_diff = c(1, -1))))

        cur_avg = cur_counts %>%
            group_by(Protocol) %>%
            summarise(average = mean(Count, na.rm = TRUE))

        cur_avg = as.data.frame(cur_avg)
        
        stat = data.frame(Contrast=paste(cur_protocols, collapse="_vs_"), CellType=celltype, 
            T=sum_cont$t.ratio, P=sum_cont$p.value, Avg1=cur_avg$average[1], Avg2=cur_avg$average[2])

        return(stat)
        
    })

    stats = as.data.frame(do.call(rbind, statList))
    stats$fdr = p.adjust(stats$P, method = "BH")

    return(stats)
    
})
## [1] 1
## [1] "Bergmann glia"
## [1] "Deep cerebellar nuclei (excitatory/inhibitory)"
## [1] "Endothelial"
## [1] "Golgi cell"
## [1] "Granule cell / lineage"
## [1] "Microglia"
## [1] "Molecular layer interneuron (basket/stellate)"
## [1] "Oligodendrocyte / OPC"
## [1] "Pericyte / SMC"
## [1] "Purkinje neuron"
## [1] "Unipolar brush cell (UBC)"
## [1] 2
## [1] "Bergmann glia"
## [1] "Deep cerebellar nuclei (excitatory/inhibitory)"
## [1] "Endothelial"
## [1] "Golgi cell"
## [1] "Granule cell / lineage"
## [1] "Microglia"
## [1] "Molecular layer interneuron (basket/stellate)"
## [1] "Oligodendrocyte / OPC"
## [1] "Pericyte / SMC"
## [1] "Purkinje neuron"
## [1] "Unipolar brush cell (UBC)"

stats = as.data.frame(do.call(rbind, statList))
stats = stats[, c("Contrast", "CellType", "T", "P", "fdr", "Avg1", "Avg2")]

stats = stats[order(stats$Contrast, stats$P), ]
stats$P = formatC(stats$P, format = "e", digits = 2)
stats$fdr = formatC(stats$fdr, format = "e", digits = 2)

datatable(
    stats,
    extensions = 'Buttons',
    caption = "Geneset Score for Cell Types",
    rownames = FALSE,
    options = list(
        dom = 'Bfrtip', # B for buttons, frtip for other table features
        buttons = c('csv', 'excel'), # Add CSV and Excel download buttons
        pageLength = 15
    )    
) %>% 
    formatRound(columns = c("T", "Avg1", "Avg2"), digits = 2)
```


```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-109314da358e6cccc5ea" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-109314da358e6cccc5ea">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"caption":"<caption>Geneset Score for Cell Types<\/caption>","data":[["SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974"],["Purkinje neuron","Bergmann glia","Unipolar brush cell (UBC)","Golgi cell","Oligodendrocyte / OPC","Pericyte / SMC","Granule cell / lineage","Endothelial","Molecular layer interneuron (basket/stellate)","Microglia","Deep cerebellar nuclei (excitatory/inhibitory)","Pericyte / SMC","Granule cell / lineage","Molecular layer interneuron (basket/stellate)","Purkinje neuron","Endothelial","Bergmann glia","Deep cerebellar nuclei (excitatory/inhibitory)","Unipolar brush cell (UBC)","Golgi cell","Oligodendrocyte / OPC","Microglia"],[-6.380269534677488,5.258810779575024,-3.810314123573261,-2.781118242052646,2.543682398018991,-2.273846788323691,-2.127126262111394,2.103183453084381,-1.019809319505402,0.7541293497893069,0.5161153052648273,-9.381343167501099,-4.745344772551483,4.243545313059662,-3.909739075581776,-4.093426217053317,3.539145651851967,3.628058971342901,3.437156229260643,2.572361890874057,2.132097461777168,1.437253637436276],["2.34e-09","6.41e-07","4.37e-04","9.00e-03","1.30e-02","2.98e-02","3.55e-02","4.34e-02","3.11e-01","4.56e-01","6.08e-01","1.06e-10","5.80e-06","5.48e-05","1.43e-04","2.69e-04","5.72e-04","6.34e-04","1.32e-03","1.49e-02","3.62e-02","1.60e-01"],["2.57e-08","3.53e-06","1.60e-03","2.48e-02","2.86e-02","5.47e-02","5.57e-02","5.97e-02","3.80e-01","5.02e-01","6.08e-01","1.16e-09","3.19e-05","2.01e-04","3.92e-04","5.92e-04","9.97e-04","9.97e-04","1.81e-03","1.83e-02","3.99e-02","1.60e-01"],[565.060825411013,3791.131949843827,656.9278663650863,1603.519728860264,902.3028894885406,193.2913262373348,726.1429321889168,106.0128697184927,1211.452881812452,18.26407247137354,1943.642981397125,193.2913262373348,726.1429321889168,1211.452881812452,565.060825411013,106.0128697184927,3791.131949843827,1943.642981397125,656.9278663650863,1603.519728860264,902.3028894885406,18.26407247137354],[4605.823570025898,929.384042348509,1437.296825996558,2569.707477623045,68.0125816990426,336.8756713909118,1174.857170536211,26.69132009879297,1411.111495135766,9.668925160283536,1761.389190936047,1140.16992884606,3839.041336083596,359.1978917469644,1294.919925861094,2804.411919253885,2048.597904183681,681.5341328512498,283.061291574286,550.4872584809983,220.3416125476707,1.553339245007982]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Contrast<\/th>\n      <th>CellType<\/th>\n      <th>T<\/th>\n      <th>P<\/th>\n      <th>fdr<\/th>\n      <th>Avg1<\/th>\n      <th>Avg2<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Bfrtip","buttons":["csv","excel"],"pageLength":15,"columnDefs":[{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":5,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":6,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[2,5,6]},{"name":"Contrast","targets":0},{"name":"CellType","targets":1},{"name":"T","targets":2},{"name":"P","targets":3},{"name":"fdr","targets":4},{"name":"Avg1","targets":5},{"name":"Avg2","targets":6}],"order":[],"autoWidth":false,"orderClasses":false,"lengthMenu":[10,15,25,50,100]}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render"],"jsHooks":[]}</script>
```


# Geneset stats for brain regions

Forebrain  
Midbrain  
Cerebellum  


``` r

regions = c("Forebrain", "Midbrain", "Cerebellum")

test_gene_data = gene_data[gene_data$Region %in% regions, ]

test_count_data = norm_data[norm_data$Gene %in% test_gene_data$MarkerGene, ]
dim(test_count_data)
## [1] 1620    7

contrasts = list()
contrasts[[1]] = c("SCTLday60", "GSE233574")
contrasts[[2]] = c("SCTLday60", "GSE247974")

nconts = length(contrasts)

statList = lapply(1:nconts, function(icont) {

    statList = lapply(regions, function(region) {

        print(region)
    
        geneset = test_gene_data[test_gene_data$Region==region, "MarkerGene"]
    
        cur_protocols = as.vector(contrasts[[icont]])
        
        cur_counts = test_count_data[which( (test_count_data$Gene %in% geneset) & (test_count_data$Protocol %in% cur_protocols) ), ]
        cur_counts$Protocol = factor(cur_counts$Protocol, levels=cur_protocols)
        cur_counts$Gene = as.factor(cur_counts$Gene)

        model = lmer(Count ~ Protocol + (1 | Gene), data = cur_counts)

        emm = emmeans(model, ~ Protocol) 
        sum_cont = summary(contrast(emm, list(sample_diff = c(1, -1))))
        
        cur_avg = cur_counts %>%
            group_by(Protocol) %>%
            summarise(average = mean(Count, na.rm = TRUE))

        cur_avg = as.data.frame(cur_avg)

        stat = data.frame(Contrast=paste(cur_protocols, collapse="_vs_"), Region=region,  
            T=sum_cont$t.ratio, P=sum_cont$p.value, Avg1=cur_avg$average[1], Avg2=cur_avg$average[2])

        return(stat)
        
    })
    
    stats = as.data.frame(do.call(rbind, statList))
    stats$fdr = p.adjust(stats$P, method = "BH")
    
    return(stats)

})
## [1] "Forebrain"
## [1] "Midbrain"
## [1] "Cerebellum"
## [1] "Forebrain"
## [1] "Midbrain"
## [1] "Cerebellum"

stats = as.data.frame(do.call(rbind, statList))
stats = stats[, c("Contrast", "Region", "T", "P", "fdr", "Avg1", "Avg2")]

stats = stats[order(stats$Contrast, stats$P), ]
stats$P = formatC(stats$P, format = "e", digits = 2)
stats$fdr = formatC(stats$fdr, format = "e", digits = 2)

datatable(
    stats,
    extensions = 'Buttons',
    caption = "Geneset Score for Region",
    rownames = FALSE,
    options = list(
        dom = 'Bfrtip', # B for buttons, frtip for other table features
        buttons = c('csv', 'excel'), # Add CSV and Excel download buttons
        pageLength = 15
    )    
) %>% 
    formatRound(columns = c("T", "Avg1", "Avg2"), digits = 2)
```


```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-6130662943320a82a602" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6130662943320a82a602">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"caption":"<caption>Geneset Score for Region<\/caption>","data":[["SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE233574","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974","SCTLday60_vs_GSE247974"],["Forebrain","Cerebellum","Midbrain","Forebrain","Cerebellum","Midbrain"],[-2.638099157429163,-1.729247273183346,-1.640168986420951,-3.227945105523987,-1.728620874639109,1.282477636796449],["8.68e-03","8.42e-02","1.02e-01","1.35e-03","8.43e-02","2.01e-01"],["2.60e-02","1.02e-01","1.02e-01","4.06e-03","1.26e-01","2.01e-01"],[655.356009654729,1234.257366008471,650.5855884414117,655.356009654729,1234.257366008471,650.5855884414117],[1015.984361000672,1551.131939474625,841.2181777218731,8506.847472619927,1516.555460336197,504.0500004492918]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Contrast<\/th>\n      <th>Region<\/th>\n      <th>T<\/th>\n      <th>P<\/th>\n      <th>fdr<\/th>\n      <th>Avg1<\/th>\n      <th>Avg2<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Bfrtip","buttons":["csv","excel"],"pageLength":15,"columnDefs":[{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":5,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":6,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[2,5,6]},{"name":"Contrast","targets":0},{"name":"Region","targets":1},{"name":"T","targets":2},{"name":"P","targets":3},{"name":"fdr","targets":4},{"name":"Avg1","targets":5},{"name":"Avg2","targets":6}],"order":[],"autoWidth":false,"orderClasses":false,"lengthMenu":[10,15,25,50,100]}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render"],"jsHooks":[]}</script>
```


# Session Information


``` r

sessionInfo()
## R version 4.4.3 (2025-02-28)
## Platform: x86_64-pc-linux-gnu
## Running under: Rocky Linux 8.7 (Green Obsidian)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/local/intel/2022.1.2.146/mkl/2022.0.2/lib/intel64/libmkl_rt.so.2;  LAPACK version 3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: America/New_York
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] emmeans_1.11.0  lme4_1.1-37     Matrix_1.7-2    car_3.1-3      
##  [5] carData_3.0-5   DT_0.33         pheatmap_1.0.12 knitr_1.50     
##  [9] ggplot2_3.5.1   dplyr_1.1.4     tidyr_1.3.1    
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.6       xfun_0.52          bslib_0.9.0        htmlwidgets_1.6.4 
##  [5] lattice_0.22-6     crosstalk_1.2.1    vctrs_0.6.5        tools_4.4.3       
##  [9] Rdpack_2.6.3       generics_0.1.3     parallel_4.4.3     pbkrtest_0.5.3    
## [13] sandwich_3.1-1     tibble_3.2.1       pkgconfig_2.0.3    RColorBrewer_1.1-3
## [17] lifecycle_1.0.4    compiler_4.4.3     munsell_0.5.1      codetools_0.2-20  
## [21] htmltools_0.5.8.1  sass_0.4.9         yaml_2.3.10        Formula_1.2-5     
## [25] pillar_1.10.1      nloptr_2.2.1       jquerylib_0.1.4    MASS_7.3-65       
## [29] cachem_1.1.0       reformulas_0.4.0   boot_1.3-31        abind_1.4-8       
## [33] multcomp_1.4-28    nlme_3.1-167       tidyselect_1.2.1   digest_0.6.37     
## [37] mvtnorm_1.3-3      purrr_1.0.4        splines_4.4.3      fastmap_1.2.0     
## [41] grid_4.4.3         colorspace_2.1-1   cli_3.6.4          magrittr_2.0.3    
## [45] survival_3.8-3     broom_1.0.8        TH.data_1.1-3      withr_3.0.2       
## [49] scales_1.3.0       backports_1.5.0    estimability_1.5.1 rmarkdown_2.29    
## [53] zoo_1.8-13         coda_0.19-4.1      evaluate_1.0.3     rbibutils_2.3     
## [57] rlang_1.1.5        Rcpp_1.0.14        xtable_1.8-4       glue_1.8.0        
## [61] rstudioapi_0.17.1  minqa_1.2.8        jsonlite_2.0.0     R6_2.6.1

rm(list=ls())

```

