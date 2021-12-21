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
png("Interaction map C42B.png")
pheatmap(log2(1 + cancer1.10kb[1:500,2:501]), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum), main = "Heatmap for C42B cell (cancer)")
dev.off()
png("Interaction map 22Rv1.png")
pheatmap(log2(1 + cancer2.10kb[1:500,2:501]), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum), main = "Heatmap for 22Rv1 cell(cancer)")
dev.off()
png("Interaction map RWPE1.png")
pheatmap(log2(1 + normal.10kb[1:500,2:501]) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum), color = colorRampPalette(c("white", "orange", "red"))(maximum), main = "Heatmap for RWPE1 cell (normal)")
dev.off()
#Subset for the regions 127Mb-131Mb.
x <- which(cancer1.10kb$V1 == "chr12:127000000:127010000")
y <- which(cancer1.10kb$V1 == "chr12:130990000:131000000")
w <- which(cancer1.40kb$V1 == "chr12:127000000:127040000")
z <- which(cancer1.40kb$V1 == "chr12:130960000:131000000")
cancer1.sub <- cancer1.10kb[x:y,(x+1):(y+1)]
cancer1.sub.40kb <- cancer1.40kb[w:z,(w+1):(z+1)]
cancer1.names.40kb <- colnames(cancer1.sub.40kb)
rm(cancer1.10kb)
rm(cancer1.40kb)

x <- which(cancer2.10kb$V1 == "chr12:127000000:127010000")
y <- which(cancer2.10kb$V1 == "chr12:130990000:131000000")
w <- which(cancer2.40kb$V1 == "chr12:127000000:127040000")
z <- which(cancer2.40kb$V1 == "chr12:130960000:131000000")
cancer2.sub <- cancer2.10kb[x:y,(x+1):(y+1)]
cancer2.sub.40kb <- cancer2.40kb[w:z,(w+1):(z+1)]
cancer2.names.40kb <- colnames(cancer2.sub.40kb)
rm(cancer2.10kb)
rm(cancer2.40kb)

x <- which(normal.10kb$V1 == "chr12:127000000:127010000")
y <- which(normal.10kb$V1 == "chr12:130990000:131000000")
w <- which(normal.40kb$V1 == "chr12:127000000:127040000")
z <- which(normal.40kb$V1 == "chr12:130960000:131000000")
normal.sub <- normal.10kb[x:y,(x+1):(y+1)]
normal.sub.40kb <- normal.40kb[w:z,(w+1):(z+1)]
normal.names.40kb <- colnames(normal.sub.40kb)
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
cancer1.40kb.norm <- vanilla_coverage(cancer1.sub.40kb)
cancer1.40kb.norm.all <- vanilla_coverage(cancer1.40kb[,2:ncol(cancer1.40kb)])
cancer2.40kb.norm <- vanilla_coverage(cancer2.sub.40kb)
cancer2.40kb.norm.all <- vanilla_coverage(cancer2.40kb[,2:ncol(cancer2.40kb)])
normal.40kb.norm  <- vanilla_coverage(normal.sub.40kb)
normal.40kb.norm.all <- vanilla_coverage(normal.40kb[,2:ncol(normal.40kb)])

#Pheatmaps for RWPE1 cell type at 40kb resolution.
maximum2 <- ceiling(max(log2(1 + normal.40kb[1:500,2:501]),log2(1 + normal.40kb.norm.all[1:500,1:500])))
png("Heatmap of an initial 40kb resolution interaction matrix (RWPE1).png")
pheatmap(log2(1 + normal.40kb[1:500,2:501]) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum2), color = colorRampPalette(c("white", "orange", "red"))(maximum2), main = "Heatmap of an initial 40kb resolution interaction matrix (RWPE1)")
dev.off()
png("Heatmap of a normalized 40kb resolution interaction matrix (RWPE1).png")
pheatmap(log2(1 + normal.40kb.norm.all[1:500,1:500]) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum2), color = colorRampPalette(c("white", "orange", "red"))(maximum2), main = "Heatmap of a normalized 40kb resolution interaction matrix (RWPE1)")
dev.off()
rm(normal.40kb)
#Bonus.

#5. Data visualization and exploration.
maximum3 <- ceiling(max(log2(1+cancer1.10kb.norm),log2(1+cancer2.10kb.norm),log2(1+normal.10kb.norm)))
png("Heatmap of a normalized C42B matrix.png")
pheatmap(log2(1 + cancer1.10kb.norm), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum3), color = colorRampPalette(c("white", "orange", "red"))(maximum3))
dev.off()
png("Heatmap of a normalized 22Rv1 matrix.png")
pheatmap(log2(1 + cancer2.10kb.norm), cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum3), color = colorRampPalette(c("white", "orange", "red"))(maximum3))
dev.off()
png("Heatmap of a normalized RWPE1 matrix.png")
pheatmap(log2(1 + normal.10kb.norm) , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(0:maximum3), color = colorRampPalette(c("white" ,"orange", "red"))(maximum3))
dev.off()
#5.3. Subtract pair of matrices.
normal.cancer1 <- sign(normal.10kb.norm-cancer1.10kb.norm)*log2(1+abs(normal.10kb.norm - cancer1.10kb.norm))
normal.cancer2 <- sign(normal.10kb.norm-cancer2.10kb.norm)*log2(1+abs(normal.10kb.norm - cancer2.10kb.norm))
cancer2.cancer1<- sign(cancer2.10kb.norm-cancer1.10kb.norm)*log2(1+abs(cancer2.10kb.norm- cancer1.10kb.norm))
maximum4 <- ceiling(max(normal.cancer1,normal.cancer2,cancer2.cancer1))
minimum  <- floor(min(normal.cancer1,normal.cancer2,cancer2.cancer1))
pheatmap(normal.cancer1 , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(minimum:maximum4), color = colorRampPalette(magma(256))(maximum4-minimum))
pheatmap(normal.cancer2 , cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(minimum:maximum4), color = colorRampPalette(magma(256))(maximum4-minimum))
pheatmap(cancer2.cancer1, cluster_rows = F, cluster_cols = F ,labels_row = '', labels_col = '', breaks = c(minimum:maximum4), color = colorRampPalette(magma(256))(maximum4-minimum))

#6. Calling of Topologically Associated Domains (TADs)
directionality_index <- function(x, w = 2e6, r = 4e4){
  bins.number <- w/r
  A <- c()
  for(i in 1:ncol(x)){
    a.range <- i-bins.number
    b.range <- i+bins.number
    if(a.range <= 0){a.range <- 1}
    if(b.range > ncol(x)){b.range <- ncol(x)}
    if(i == 1){
      a <- 0
      b <- sum(x[i,(i+1):b.range])
      e <- (a+b)/2
      DI <-(((b-a)/abs(a-b))*((((a-e)^2)/e) + (((b-e)^2)/e)))
    }else if(i == ncol(x)){
      a <- sum(x[a.range:(i-1),i])
      b <- 0
      e <- (a+b)/2
      DI <-(((b-a)/abs(a-b))*((((a-e)^2)/e) + (((b-e)^2)/e)))
    }else{
      a <- sum(x[a.range:(i-1),i])
      b <- sum(x[i,(i+1):b.range])
      e <- (a+b)/2
      DI <-(((b-a)/abs(a-b))*((((a-e)^2)/e) + (((b-e)^2)/e)))
    }
    A[i] <- DI
  }
  return(A)
}

cancer1.di<- directionality_index(cancer1.40kb.norm.all)
cancer2.di<- directionality_index(cancer2.40kb.norm.all)
normal.di <- directionality_index(normal.40kb.norm.all)
cancer1.df <- data.frame(position = cancer1.40kb$V1 ,DI = cancer1.di)
cancer2.df <- data.frame(position = cancer2.40kb$V1 ,DI = cancer2.di)
normal.df  <- data.frame(position = normal.40kb$V1 ,DI = normal.di)

w <- which(cancer1.df$position == "chr12:127000000:127040000")
z <- which(cancer1.df$position == "chr12:130960000:131000000")
cancer1.sub.df <- cancer1.df[w:z,2]

w <- which(cancer2.df$position == "chr12:127000000:127040000")
z <- which(cancer2.df$position == "chr12:130960000:131000000")
cancer2.sub.df <- cancer2.df[w:z,2]

w <- which(normal.df$position == "chr12:127000000:127040000")
z <- which(normal.df$position == "chr12:130960000:131000000")
normal.sub.df <- normal.df[w:z,2]

xmin.cancer1 <- c()
xmax.cancer1 <- c()
for (name in strsplit(cancer1.names.40kb,":")){
  xmin.cancer1 <- c(xmin.cancer1,strtoi(name[2]))
  xmax.cancer1 <- c(xmax.cancer1,strtoi(name[3]))
}

xmin.cancer2 <- c()
xmax.cancer2 <- c()
for (name in strsplit(cancer2.names.40kb,":")){
  xmin.cancer2 <- c(xmin.cancer2,strtoi(name[2]))
  xmax.cancer2 <- c(xmax.cancer2,strtoi(name[3]))
}

xmin.normal <- c()
xmax.normal <- c()
for (name in strsplit(normal.names.40kb,":")){
  xmin.normal <- c(xmin.normal,strtoi(name[2]))
  xmax.normal <- c(xmax.normal,strtoi(name[3]))
}
df.1 <- data.frame(xmin.cancer1 = xmin.cancer1, xmax.cancer1= xmax.cancer1, ymin=0,ymax.1 = cancer1.sub.df, negative = as.factor(sign(cancer1.sub.df)))
df.2 <- data.frame(xmin.cancer2 = xmin.cancer2, xmax.cancer2= xmax.cancer2, ymin=0,ymax.2 = cancer2.sub.df, negative = as.factor(sign(cancer2.sub.df)))
df.n <- data.frame(xmin.normal  = xmin.normal , xmax.normal = xmax.normal , ymin=0,ymax.n = normal.sub.df , negative = as.factor(sign(normal.sub.df)))
png("DI for C42B cell type.png")
ggplot(df.1) + geom_rect(aes(xmin = xmin.cancer1 ,xmax = xmax.cancer1,ymin = ymin,ymax = ymax.1, col=negative)) + ggtitle("DI for C42B cell type") + xlab("Chromosome position") + ylab("DI") + theme(legend.position = "None")
dev.off()
png("DI for 22Rv1 cell type.png")
ggplot(df.2) + geom_rect(aes(xmin = xmin.cancer2 ,xmax = xmax.cancer2,ymin = ymin,ymax = ymax.2, col=negative)) + ggtitle("DI for 22Rv1 cell type")+ xlab("Chromosome position") + ylab("DI") + theme(legend.position = "None")
dev.off()
png("DI for RWPE1 cell type.png")
ggplot(df.n) + geom_rect(aes(xmin = xmin.normal  ,xmax = xmax.normal ,ymin = ymin,ymax = ymax.n, col=negative)) + ggtitle("DI for RWPE1 cell type")+ xlab("Chromosome position") + ylab("DI") + theme(legend.position = "None")
dev.off()
