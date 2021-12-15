#Set the environment.
library("ggplot2")
library("data.table")
library("pheatmap")
library("viridisLite")
library("viridis")
setwd("C:/Users/Daniel/OneDrive/Bureau/EPFL/Génétique & Génomique/Genetics_project_2")
#Load the data.
#10kb resolution.
cancer1.10kb = fread("C42B_chr12_10kb_hic_matrix.txt")
cancer2.10kb = fread("22Rv1_chr12_10kb_hic_matrix.txt")
normal.10kb  = fread("RWPE1_chr12_10kb_hic_matrix.txt")
#40kb resolution.
cancer1.40kb = fread("C42B_chr12_40kb_hic_matrix.txt")
cancer2.40kb = fread("22Rv1_chr12_40kb_hic_matrix.txt")
normal.40kb  = fread("RWPE1_chr12_40kb_hic_matrix.txt")

#Pheatmaps.
maximum <- ceiling(max(log2(1 + cancer1.10kb[1:500,2:501]),log2(1 + cancer2.10kb[1:500,2:501]),log2(1 + normal.10kb[1:500,2:501])))
pheatmap(log2(1 + cancer1.10kb[1:500,2:501]), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum))
pheatmap(log2(1 + cancer2.10kb[1:500,2:501]), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum))
pheatmap(log2(1 + normal.10kb[1:500,2:501]) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum))

#Subset for the regions 127Mb-131Mb.
x <- which(cancer1.10kb$V1 == "chr12:127000000:127010000")
y <- which(cancer1.10kb$V1 == "chr12:130990000:131000000")
cancer1.sub <- cancer1.10kb[x:y,(x+1):(y+1)]
rm(cancer1.10kb)

x <- which(cancer2.10kb$V1 == "chr12:127000000:127010000")
y <- which(cancer2.10kb$V1 == "chr12:130990000:131000000")
cancer2.sub <- cancer2.10kb[x:y,(x+1):(y+1)]
rm(cancer2.10kb)

x <- which(normal.10kb$V1 == "chr12:127000000:127010000")
y <- which(normal.10kb$V1 == "chr12:130990000:131000000")
normal.sub <- normal.10kb[x:y,(x+1):(y+1)]
rm(normal.10kb)

#Implement the Vanilla Coverage normalization.
vanilla_coverage <- function(x){
  sum.the.rows <- rowSums(x)
  by.zero <- which(sum.the.rows == 0)
  if(length(by.zero) != 0){
    for(i in by.zero){
      specific.norm <- diag(1/sum.the.rows)
      specific.norm[i,i] = 0
    }
  }else{specific.norm <- diag(1/sum.the.rows)}
  norm.matrix <- (specific.norm %*% as.matrix(x) %*% specific.norm)
  total.sum <- sum(rowSums(x))
  norm.total.sum <- sum(rowSums(norm.matrix))
  return(norm.matrix*(total.sum/norm.total.sum))
}

#Normalization
cancer1.10kb.norm <- vanilla_coverage(cancer1.sub)
cancer2.10kb.norm <- vanilla_coverage(cancer2.sub)
normal.10kb.norm  <- vanilla_coverage(normal.sub)
cancer1.40kb.norm <- vanilla_coverage(cancer1.40kb[,2:ncol(cancer1.40kb)])
rm(cancer1.40kb)
cancer2.40kb.norm <- vanilla_coverage(cancer2.40kb[,2:ncol(cancer2.40kb)])
rm(cancer2.40kb)
normal.40kb.norm  <- vanilla_coverage(normal.40kb[,2:ncol(normal.40kb)])

#Pheatmaps for RWPE1 cell type at 40kb resolution.
maximum2 <- ceiling(max(log2(1 + normal.40kb[1:500,2:501]),log2(1 + normal.40kb.norm[1:500,1:500])))
pheatmap(log2(1 + normal.40kb[1:500,2:501]) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum2), color = colorRampPalette(c("white", "orange", "red"))(maximum2))
pheatmap(log2(1 + normal.40kb.norm[1:500,1:500]) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum2), color = colorRampPalette(c("white", "orange", "red"))(maximum2))
rm(normal.40kb)
#Bonus.

#5. Data visualization and exploration.
maximum3 <- ceiling(max(log2(1+cancer1.10kb.norm),log2(1+cancer2.10kb.norm),log2(1+normal.10kb.norm)))
pheatmap(log2(1 + cancer1.10kb.norm), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum3), color = colorRampPalette(c("white", "orange", "red"))(maximum3))
pheatmap(log2(1 + cancer2.10kb.norm), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum3), color = colorRampPalette(c("white", "orange", "red"))(maximum3))
pheatmap(log2(1 + normal.10kb.norm) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum3), color = colorRampPalette(c("white" ,"orange", "red"))(maximum3))

#5.3. Subtract pair of matrices.
normal.cancer1 <- normal.10kb.norm - cancer1.10kb.norm
normal.cancer2 <- normal.10kb.norm - cancer2.10kb.norm
cancer1.cancer2<- cancer1.10kb.norm- cancer2.10kb.norm
maximum4 <- ceiling(max(1+normal.cancer1,1+normal.cancer2,1+cancer1.cancer2))
minimum  <- floor(min(1+normal.cancer1,1+normal.cancer2,1+cancer1.cancer2))
pheatmap(normal.cancer1 , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(minimum:maximum4), color = colorRampPalette(magma(256))(maximum4-minimum))
pheatmap(normal.cancer2 , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(minimum:maximum4), color = colorRampPalette(magma(256))(maximum4-minimum))
pheatmap(cancer1.cancer2, cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(minimum:maximum4), color = colorRampPalette(magma(256))(maximum4-minimum))

#6. Calling of Topologically Associated Domains (TADs)
directionality_index <- function(x){
  
  e <- (a+b)/2
  DI <- (((a+b)/abs(a-b))*(((a-e)^2)/e + ((b-e)^2)/e))
  return(DI)
}