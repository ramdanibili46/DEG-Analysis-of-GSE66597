#Modul: Analisis Infeksi Virus H5N1 terhadap Ekspresi Gen Astrosi Manusia
#Dataset: GSE 66597
#Platform: GPL6480
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)

# Memanggil library
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# Pengambilan data dari GEO
#GEO (Gene Expression Omnibus) adalah database publik milik NCBI
gset <- getGEO("GSE66597", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Membuat nama kolom yang tepat agar sesuai dengan variabel top tabel
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Pengelompokkan berdasarkan jenis sampel
gsms <- "000000000111111111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
# log2 digunakan untuk:
#1. Menstabilkan varians
#2. Mendekati asumsi model linear
#3. Memudahkan interpretasi log fold change
# log transform adalah variabel logika
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# Memasukkan kelompok sampel dan mendefinisikan kelompok sampel
gs <- factor(sml)
groups <- make.names(c("Control","H5N1"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest dan mengkalkulasi kembali model koefisien
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualisasi dan quality control hasil tes
# Membangun histogram dari P-values untuk semua genes. Normal test
# Asumsi tidak semua gen merupakan ekspresi diferensial.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")

# Merangkum hasil sebagai "up", "down" atau "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# Membuat Q-Q plot untuk t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest

# volcano plot menggunakan limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# Analisis ekspresi data umum
ex <- exprs(gset)

# Visualisasi Data
# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE66597", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE66597", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 8, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=8", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), 
       legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", 
       pt.cex=1.5)