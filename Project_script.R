library("ggplot2")
library("data.table")
library("pheatmap")
setwd("C:/Users/Daniel/OneDrive/Bureau/EPFL/Génétique & Génomique/Genetics_project_2")

cancer1 = fread("C42B_chr12_10kb_hic_matrix.txt")
cancer2 = fread("22Rv1_chr12_10kb_hic_matrix.txt")
normal = fread("RWPE1_chr12_10kb_hic_matrix.txt")

maximum <- max(log2(1 + cancer1[1:500,2:501]),log2(1 + cancer2[1:500,2:501]),log2(1 + normal[1:500,2:501]))

pheatmap(log2(1 + cancer1[1:500,2:501]), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum))
pheatmap(log2(1 + cancer2[1:500,2:501]), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum))
pheatmap(log2(1 + normal[1:500,2:501]), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum))


cancer1.sub <- cancer1[]
specific.norm <- diag(1/rowSums(cancer1["chr12:127000:128000":"chr12:130000:131000",2:501]))
matr <- (specific.norm %*% as.numeric(cancer1[1:500,2:501]) %*% specific.norm)
