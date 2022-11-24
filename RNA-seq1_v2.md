---
title: "RNA-seq1"
author: Siwen Chen, Midhuna Immaculate Joseph Maran, Sanjay Kumar Srikakulam
date:
theme: united
highlight: tango
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
    toc_float:
      collapsed: true
      smooth_scroll: false
tags: rnaseq
---


# Get to know the system

1. How many nodes, and how many CPUs and how much RAM does each node have
    - 4 nodes
    - 12 CPUs and 121000 MiB RAM per node
2. What is the volume’s size of `/vol/COMPEPIWS`
    - 3.9 TiB
3. Temp folder location
    - `/vol/COMPEPIWS/tmp`
4. Group working directory
    - `/vol/COMPEPIWS/groups/rnaseq1`
5. 
6. FASTQ files (data)
    - `ls /vol/COMPEPIWS/data/reduced/RNA-seq/`
    - 



# Tips and Tricks
1. `ssh -X -A -i ~/.ssh/id_rsa seven_epigen@193.175.249.28`
2. `ssh -A bibigrid-worker-1-3-etqcgckjrdsvsxe`
3. `screen -s screen_name`
4. `Ctrl-a c`
5. `Ctrl-a c`
6. `Ctrl-a d`
7. `Ctrl-d`
8. `Ctrl-d`
9. `ssh -X -A -i ~/.ssh/id_rsa seven_epigen@193.175.249.28`
10. `ssh -A bibigrid-worker-1-3-etqcgckjrdsvsxe`
11. `screen -ls`
12. `screen -rd screen_name`
13. `Ctrl-a n` and `Ctrl-a p`


# Working with tables on the command line, or awk

3. `'wc -l /vol/COMPEPIWS/pipelines/references/genome_genes.gtf'`
    - 671462 lines

4. `'grep -c "exon" /vol/COMPEPIWS/pipelines/references/genome_genes.gtf'`
    - 324748 exons

    or `'awk '/exon/{++cnt} END {print cnt}' /vol/COMPEPIWS/pipelines/references/genome_genes.gtf'`
    - 324748 exons

5. `'awk '(/exon/ && $5-$4 > 1000) {++cnt} END {print cnt}' /vol/COMPEPIWS/pipelines/references/genome_genes.gtf'`
    - 19593

6. `'awk '(/exon/ && $10 ~ /Sox17/) {++cnt} END {print cnt}' /vol/COMPEPIWS/pipelines/references/genome_genes.gtf'`
    - 20

7. `'awk '(/exon/ && $1 ~ /chr2/) {++cnt} END {print cnt}' /vol/COMPEPIWS/pipelines/references/genome_genes.gtf'`
    - 29373

8. `'awk -F" " '{print $3}' /vol/COMPEPIWS/pipelines/references/genome_genes.gtf | sort | uniq -c'`
    - 286783 CDS
    - 324748 exon
    - 29978 start_codon
    - 29953 stop_codon

    or `'awk '{count[$3]++} END {for (word in count) print word, count[word]}'  /vol/COMPEPIWS/pipelines/references/genome_genes.gtf'`
    - exon 324748
    - CDS 286783
    - start_codon 29978
    - stop_codon 29953

9. `'sort -t " " -k1 -k3  /vol/COMPEPIWS/pipelines/references/genome_genes.gtf >> sorted_genomes.gtf'`

    or `'sort -k 1 -k 3 /vol/COMPEPIWS/pipelines/references/genome_genes.gtf > genome_genes_new_file.gtf'`


# Conda and Bioconda

a. `source /vol/COMPEPIWS/conda/miniconda3/bin/activate`
b. `conda create -p /vol/COMPEPIWS/groups/rnaseq1/conda/<srikakulam>_rnaseq1`
c. `conda activate /vol/COMPEPIWS/groups/rnaseq1/conda/<srikakulam>_rnaseq1`
d. `conda install -c r r-ggplot2`, `conda install -c bioconda fastqc bedtools==2.22`
f. `conda update -c bioconda bedtools`, version `conda list | grep bedtools` = 2.30.0
g. `conda deactivate`

# Basics in R 
## Data frames
``` {r}
# 2. Load required libraries
library(ggplot2)
library(reshape2)
library(GenomicRanges)

# 3. Load tables
kidney = read.table("/vol/COMPEPIWS/groups/shared/kidney_14.5_mouse_1_2kbW_bed_counts.txt", header=TRUE)
liver = read.table("/vol/COMPEPIWS/groups/shared/liver_14.5_mouse_1_2kbW_bed_counts.txt", header=TRUE)

# 4. Dimension of each table
dim(kidney)     # => 1362728       6
dim(liver)      # => 1362728       6

# 5. What are the column names of each table
colnames(kidney)    # => [1] "Chr"      "Start"    "End"      "H3K27me3" "H3K36me3" "H3K9me3"
colnames(liver)     # => [1] "Chr"      "Start"    "End"      "H3K27me3" "H3K36me3" "H3K9me3"

# 6. Calculate the genome length from each data set
liver[nrow(liver), "End"] - liver[1, "Start"]       # => 91742000
kidney[nrow(kidney), "End"] - kidney[1, "Start"]    # => 91742000

# 7. Create a data frame by concatenating vertically the two data sets. Add a new column called (cell_type) to annotate the rows of liver table with "liver" and the kidney’s rows with "kidney"
kidney["cell_type"] = "kidney"
liver["cell_type"] = "liver"
cell_data = rbind(kidney, liver)

# 8. What is the dimension of the new data frame?
dim(cell_data)      # => 2725456       7

# 9. Reshape the data frame (excluding chr, start, end) from wide into long format using "cell_type" as variable id
cell_data_melt = melt(cell_data[4:7], id=c("cell_type"))
dim(cell_data_melt)     # => 8176368       3
# 10.  
ggplot(cell_data_melt, aes(x = value, color = variable, fill=variable)) + geom_density(alpha=0.3) + facet_wrap(~cell_type)
```

 ![](https://i.imgur.com/TNmgt5p.png)

``` R
# 11.   

ggplot(cell_data_melt, aes(x = value, color = variable, fill=variable)) + geom_density(alpha=0.3) + facet_wrap(~cell_type)+xlim(0,100)
```
![](https://i.imgur.com/1Pg53VW.png)
``` R
# 12. 

hori_merge <- data.frame(liver, kidney)
horizondal <- data.frame(hori_merge$H3K27me3, hori_merge$H3K36me3, hori_merge$H3K9me3, hori_merge$H3K27me3.1, hori_merge$H3K36me3.1, hori_merge$H3K9me3.1)
# 13. 
colnames(horizondal) <- c("H3K27me3_liver","H3K36me3_liver","H3K9me3_liver","H3K27me3_kidney","H3K36me3_kidney","H3K9me3_kidney")
# 14. 
pca_hori <- prcomp(horizondal, scale=TRUE)
# 15. 
summary(pca_hori)  #0.5063 0.2328
# 16. 
png("screeplot.PNG", width=16,height=10,units="in",res=600)
    par(mar=c(7.1, 7.1, 7.1, 3.1))
    screeplot(pca_hori)
    dev.off()
    
pcaCharts <- function(x) {
    x.var <- x$sdev ^ 2
    x.pvar <- x.var/sum(x.var)
    print("proportions of variance:")
    print(x.pvar)
    
    par(mfrow=c(2,2))
    plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
    plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
    screeplot(x)
    screeplot(x,type="l")
    par(mfrow=c(1,1))
 
png("cumu_screeplot.png",width=16,height=10,units="in",res=600)
par(mar=c(7.1, 7.1, 7.1, 3.1))
pcaCharts(s)
dev.off()
```
![](https://i.imgur.com/YqmA8ks.png)
``` R
# 17.
head(pca_hori$rotation) #6

# 18.
values <- as.data.frame(pca_hori$rotation)
ggplot(data = values, aes(x = PC1, y = PC2, label = rownames(values))) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65") + geom_text(colour = "tomato", alpha = 0.8, size = 4) +ggtitle("PCA- liver and kidney")
```
![](https://i.imgur.com/guW6KzI.png)






## GenomicRanges

``` R
# 1. 
kidney_grange <- makeGRangesFromDataFrame(kidney, keep.extra.columns=TRUE)
liver_grange <- makeGRangesFromDataFrame(liver, keep.extra.columns=TRUE)

# 2. 
sum(width(reduce(kidney_grange)))
# 2725456021 bases
sum(width(reduce(liver_grange)))
# 2725456021 bases

# 3. 
sub_kidney_grange <- kidney_grange[seqnames(kidney_grange) == "chr2"]

# 4. 
up100_kidney_grange <- flank(kidney_grange, 100)

# 5. 
overlapping <- intersect(kidney_grange, liver_grange)

# 6. 
grl <- GRangesList("txA" = kidney_grange, "txB" = liver_grange)
```

# SLURM
1. `mkdir <username>_slurm_task`
2. `squeue`
3. `squeue -u <username>`
4. `for i in $(seq 1 15)
    do
        sbatch -J test_${i} -o output_${i}.out --wrap="sleep 20 && hostname"
    done`
5. `squeue`

# Nextflow RNA-seq pipeline

1. Organize your raw data
    1.b. ls /vol/COMPEPIWS/data/reduced/RNA-seq/
    1.c. 2 cell types, 2 timepoints each, 2 replicates each

2. Generate the samplesheet
    1. `cd /vol/COMPEPIWS/groups/rnaseq1/tasks`
    2. samplesheet.csv

        ``` bash
        echo "sample,fastq_1,fastq_2,strandedness" > samplesheet.csv
        for i in /vol/COMPEPIWS/data/reduced/RNA-seq/*.gz; do 
            file_name="$i"
            sample_id=$(basename $i | cut -f 1,2 -d '.')
            strandness='unstranded'
            echo "$sample_id,$file_name,,$strandness" >> samplesheet.csv
        done
        ```
3. Running the pipeline (from the master node)
    1. Activate the conda env
        - `source /vol/COMPEPIWS/conda/miniconda3/bin/activate /vol/COMPEPIWS/conda/miniconda3/envs/core/`
    2. Run the pipeline

        ```bash
        nextflow run nf-core/rnaseq -r 3.2 -profile singularity -c /vol/COMPEPIWS/pipelines/configs/rnaseq.config --input /vol/COMPEPIWS/groups/rnaseq1/tasks/samplesheet.csv --aligner star_rsem
        ```
# Quality control

## Exploring the result folder:

1. - 10m 30s, NFCORE_RNASEQ:RNASEQ:QUANTIFY_RSEM:RSEM_CALCULATEEXPRESSION (3.9 GB)
2. /vol/COMPEPIWS/groups/rnaseq1/tasks/results/star_rsem
3. find . -name "*.bigWig"
    - vol/COMPEPIWS/groups/rnaseq1/tasks/results/star_rsem/bigwig
4. /vol/COMPEPIWS/groups/rnaseq1/tasks/results/star_rsem
    - gene_id, transcript_id(s), length, effective_length, expected_count, TPM, FPKM
    - contains gene level expression estimates

## MultiQC

### FastQC(Raw):

1. No. of reads:
    - kidney_14.5_RNA_1: 2242537
    - kidney_14.5_RNA_2: 2060767
    - kidney_15.5_RNA_1: 2722856
    - kidney_15.5_RNA_2: 2969779
    - liver_14.5_RNA_1: 2156312
    - liver_14.5_RNA_2: 1584391
    - liver_15.5_RNA_1: 1926407
    - liver_15.5_RNA_2: 2530826

2. Quality of bases per reads:
    - Good, above PHRED score 30. All samples passed the QC

3. Quality of reads:
    - Good, mean around 38

4. Read lengths:
    - 20-100 bases

5. Adapter: 
    - Adapter contamination? yes 
    - Adapters were eliminated after trimming

6. After trimming:
    - Negligible number of reads were lost

### Samtools:

1. Mapping efficiency:
    - Good. Over 98% of the sequcences in each samples were mapped.

2. Duplication rate (%): 
    Over 50% duplication rate
    - kidney_14.5_RNA_1: 52
    - kidney_14.5_RNA_2: 53.4
    - kidney_15.5_RNA_1: 57.3
    - kidney_15.5_RNA_2: 57.6
    - liver_14.5_RNA_1: 67.3
    - liver_14.5_RNA_2: 63.4
    - liver_15.5_RNA_1: 62.3
    - liver_15.5_RNA_2: 66.2

### Rsem:
Many multimapping reads?
    - No, most reads are aligned within 6 reference regions, only few reads have more than 6 alignments.
    - Most of the reads were mapped to 1 reference regions; less than 500000 reads map to more than 1 reference regions.

### RSeQC:
1. Read distribution:
    - About 65 - 80% of the reads map to exonic regions in which about 50% of the reads map to CDS exonic regions. The majority of the data of exonic and only a negligible percent of reads mapping to other intergenic regions. Therefore, the data is a good library. 

2. Library strandedness: 
    - Since the ratio of sense and antisense are almost 50%:50%, the libraries are unstranded (if the ratio is rather higher or lower, the libraries will be stranded).

### DupRadar:
The duplication rate and library complexity:
    - High duplication rate (~94%) happens when the genes have high expression, and high duplication numbers at high read counts, hence the technical duplication is little and the library is with high complexity. 


# IGV

We chose gene cndp2 (Chr18) as our gene of interest and carried out organ-specific grouping and analysis of gene expression.

![](https://i.imgur.com/DkL3Cbv.png)
The image shows the gene expression for kidney samples

![](https://i.imgur.com/pQRQt9X.png)
![](https://i.imgur.com/e8Q2F45.png)

We grouped the replicates (kidney) together to compare the expression levels. The difference in expression levels among the replicates' forward and reverse strands are minimal.

![](https://i.imgur.com/Wb3ZSrE.png)
![](https://i.imgur.com/zWKpe2b.png)

We also grouped and analyzed the samples based on timepoints. Here, we see a moderately significant increase in expression among the 15.5 samples (kidney)


![](https://i.imgur.com/cYQblDX.png)

The same grouping was carried out for liver samples. Shown is the screenshot of expression levels in liver for the same gene.

![](https://i.imgur.com/Iplbltt.png)
![](https://i.imgur.com/PiTMQJw.png)
Replicates of 15.5 timepoints show a difference in expression levels (probably a biological variance?)


![](https://i.imgur.com/OuFBegK.png)
![](https://i.imgur.com/1JtV4UQ.png)

Similarly, time-point based differences are seen though not consistent with both the replicates. 

# Exploration and Differential analysis with edgeR

``` R
# load library
library(edgeR)
library(tximport)
library(biomaRt)
library(ggplot2)
library(ComplexHeatmap)
library(stringi)

#step 5 - Creating an annotation and loading data
# a)
genes <- list.files("/vol/COMPEPIWS/groups/rnaseq1/tasks/results/star_rsem/", pattern="*genes.results", all.files=FALSE, full.names=TRUE)
isoforms <- list.files("/vol/COMPEPIWS/groups/rnaseq1/tasks/results/star_rsem/", pattern="*isoforms.results", all.files=FALSE, full.names=TRUE)
sampleID = stri_sub(genes, from=56, to=-15)
timepoint <- as.numeric(stri_sub(genes, from=-24, to=-21))
replicate <- as.numeric(stri_sub(genes, from=-15, to=-15))
group <- stri_sub(genes, from=56, length=6)
group <- gsub("_", "",group)
df <- data.frame(genes, isoforms, sampleID, timepoint, replicate, group)

# b)
tmp <- read.table(genes[1], header = TRUE, sep='\t')
tx2gene <- tmp[, c("transcript_id.s.", "gene_id")]
txi <- tximport(genes, type="rsem", tx2gene=tx2gene)

# step 6 - Design matrix and DGEList objext
design_matrix = model.matrix(~group, data = txi)
dge_list = DGEList(counts=txi$counts)
nrow(dge_list)     # Number of genes is 1392

# 7. Filtering
to_keep = rowSums(cpm(dge_list) > 0.5) >= 3
filtered_dge_list = dge_list[to_keep, ,keep.lib.sizes=FALSE]
dim(filtered_dge_list)      # Number of genes after filtering is 1127, so totally 265 genes were filtered

# 8. PCA (PC1 explains 86% and PC2 explains 9% of the variance)

## update to PCA plot 

png("q8_pca.png",width=16,height=10,units="in",res=600)
par(mar=c(7.1, 7.1, 7.1, 3.1))
plotMDS(filtered_dge_list$counts, top=1127, plot=TRUE, labels=sampleID, xlim=c(-3200,3200), ylim=c(-1000, 1000), gene.selection="common")
dev.off()
```

![](https://i.imgur.com/KYyexAU.png)


``` R
# 9. Normalization: Normalize the data (uses TMM method by default)
eff_lib = calcNormFactors(filtered_dge_list)

# 10. Differential analysis:
# Quasi-likelihood F test:
y <- estimateDisp(eff_lib,design_matrix)
fit <- glmQLFit(y,design_matrix)
qlf <- glmQLFTest(fit,coef=2)
png("logCPM.PNG", width=16,height=10,units="in",res=600)
par(mar=c(7.1, 7.1, 7.1, 3.1))
plotMD(qlf)
dev.off()
```
![](https://i.imgur.com/2gStwDt.png)

``` R
# Table of the Top Differentially Expressed Genes/Tags
top_genes = topTags(qlf, n=nrow(qlf))

# Data frame with results for all genes
qlf_coeff = as.data.frame(qlf$coefficients)

# 11. BiomaRt
# a) 

#httr::set_config(httr::config(ssl_verifypeer = FALSE)) # in case of the connection problem, we can use this line here
biomartConnection<-useMart("ensembl",dataset="mmusculus_gene_ensembl")
geneAnnotation<-getBM(attributes=c("ensembl_gene_id", "external_gene_name", "description","chromosome_name", "start_position", "end_position"),mart=biomartConnection)

# list all attributes of biomartConnection
attr_biomart = listAttributes(biomartConnection)
listAttributes(biomartConnection)[1:5,] # show the first 5 attributes, overall there are 2974 attributes.

#b)
top_genes_df = as.data.frame(top_genes)
top_genes_df["gene_names"] = rownames(top_genes)
top_gene_Annotation = merge(top_genes_df, geneAnnotation, by.x="gene_names", by.y="external_gene_name") 
# loss 68 genes here, not all gene names exit in external_gene_name - 1062 genes left. 
# Also, check duplicates in top_gene_Annotation by gene_names, 3 duplicates, only their ensembl_gene_id are different. 
# the duplicates are from geneAnnotation, and excluding different ensemble id, the other values are all the same.

# 12. Create a bedgraph file

# add values in the value position
log10_fdr = -log10(top_gene_Annotation$FDR)
log10_fdr[log10_fdr == 0] <- 1e-10
top_gene_Annotation$value <- log10_fdr


# create a bedGraph file - 1062 genes
top_gene_Annotation_bedgraph = top_gene_Annotation[, 9:12]
write.table(top_gene_Annotation_bedgraph,"DEG_table.bedGraph",sep="\t",col.names=FALSE,row.names=FALSE, quote=FALSE)

# 13. CPM file - 1059 genes
cpm_val = cpm(filtered_dge_list$counts)
cpm_val = as.data.frame(cpm_val)
colnames(cpm_val) = sampleID
cpm_val["gene_names"] = rownames(cpm_val)

# excluding the non-geneAnnotation genes 
cpm_val1 = merge(cpm_val, top_gene_Annotation, by.x = "gene_names", by.y = "gene_names")
cpm_val = cpm_val1[1:9]

write.table(cpm_val,"CPM_expression.csv",sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

# 14. Exporting
meanCPMKidney = as.data.frame(apply(cpm_val[, 2:5], 1, mean))
colnames(meanCPMKidney) = "meanCPMKidney"
meanCPMKidney$gene_names = cpm_val$gene_names
meanCPMLiver = as.data.frame(apply(cpm_val[, 6:9], 1, mean))
colnames(meanCPMLiver) = "meanCPMLiver"
meanCPMLiver$gene_names = cpm_val$gene_names

# a)
bed_file = data.frame(top_gene_Annotation$chromosome_name, top_gene_Annotation$start_position, top_gene_Annotation$end_position, top_gene_Annotation$gene_names, top_gene_Annotation$value, top_gene_Annotation$logFC, meanCPMKidney$meanCPMKidney, meanCPMLiver$meanCPMLiver, top_gene_Annotation$ensembl_gene_id)
colnames(bed_file) = c("Chromosome", "Start", "End", "gene_ID", "-log10FDR", "logFC", "meanCPMKidney", "meanCPMLiver", "EnsembleID")
bed_file$Chromosome = paste("chr", bed_file$Chromosome, sep="")

# save the result
write.table(bed_file,"DEGs_complete.bed",sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

# b)
# select the significant genes
library(dplyr)
DEGs_sig = bed_file %>% filter((top_gene_Annotation$FDR <= 0.01) & abs(top_gene_Annotation$logFC >=2))

# save the result
write.table(DEGs_sig,"DEGs_sig_01FDR_2FC.bed",sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

#For IGV:
write.table(bed_file, "DEGs_complete_IGV.bed", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(DEGs_sig, "DEGs_sig_01FDR_2FC_IGV.bed", sep = "\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

# 15. Volcano plot:
top_volcano <- top_gene_Annotation
top_volcano$significant <- "NO"
top_volcano$significant[top_gene_Annotation$FDR <= 0.01 & abs(top_gene_Annotation$logFC) >= 2] = "YES"

png("Volcano_plot.PNG", width=16,height=10,units="in",res=600)
par(mar=c(7.1, 7.1, 7.1, 3.1))
ggplot(data=top_volcano,aes(x=logFC, y=value,col=significant)) +geom_point() + theme_minimal()+ geom_vline(xintercept=c(-2, 2), col="red") +
    geom_hline(yintercept=2, col="red")
dev.off()

```
![](https://i.imgur.com/WN2VBeL.png)


```R
# 16. Heatmap:

# means dataframe
means_heatmap = data.frame(DEGs_sig$gene_ID, DEGs_sig$meanCPMKidney, DEGs_sig$meanCPMLiver)
colnames(means_heatmap) = c("gene_ID", "meanCPMKidney", "meanCPMLiver")
rownames(means_heatmap) <- means_heatmap$gene_ID
means_heatmap <- means_heatmap[,-1]

# draw and save heatmap
means_heatmap <- data.matrix(means_heatmap)
png("heatmap.PNG", width=16,height=10,units="in",res=600)
par(mar=c(7.1, 7.1, 7.1, 3.1))
pheatmap(means_heatmap, fontsize=5)
dev.off()
```

![](https://i.imgur.com/mAA0GNS.png)

# Integrative analysis

## 2 Integrative data exploration in IGV

![](https://i.imgur.com/GxfhVYw.png)
![](https://i.imgur.com/NSCUuCe.png)
![](https://i.imgur.com/bDdquIJ.png)
![](https://i.imgur.com/QqXXcvy.png)
![](https://i.imgur.com/vPQgfj4.png)
![](https://i.imgur.com/jm6yGHP.png)
Expression of cndp2 more in kidney than liver; at 15.5 than 14.4 in both organ types
![](https://i.imgur.com/YZMIvZf.png)
Kidney_14.5(cndp2)
Narrow peaks: H3k4me3, H3k9ac
Broad peaks: H3k4me1, H3k36me3
![](https://i.imgur.com/A4m2PWf.png)
Kidney_15.5 hr
Narrow peaks: H3k4me3
Broad peaks:H3k4me1, H3k36me3
![](https://i.imgur.com/iS1ulM0.png)
liver_14.5 hr
Narrow peak: H3k4me3,
Broad peaks: H3k4me1, H3k27ac, H3k36me3
![](https://i.imgur.com/QKEouI3.png)
Liver_15.5 hr
Narrow peak: H3k4me3
Broad peaks:H3k4me1, H3k36me3
![](https://i.imgur.com/QA39ueP.png)
DEGene: GDA, Kidney_14.5 hr
Narrow peaks: H3k4me3, H3k9ac
Broad peaks:H3k4me1, H3k36me3, H3k27me3
Possible Polycomb repression

## 3 Integrative analysis

### 1 Summarize signals over regions of interest

``` R
# a)
library(data.table)
enhancers <- read.table("/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/segmentation/kidney_15_segments_dense.bed", head=FALSE, sep="\t", fill=TRUE)
enhancers <- as.data.frame(enhancers)
enhancers <- enhancers[-1,]
enhancers <- enhancers[enhancers$V4 %like% "Enh",]

promoters <- read.table("/vol/COMPEPIWS/pipelines/references/mm10_genome_genes.bed", head=FALSE, sep="\t")




```

### 6 Chromatin states and gene expression

``` R
# a)
# only pick the kidney_14.5_rep1 file
gene_result <- as.data.frame(read.table("/vol/COMPEPIWS/groups/rnaseq1/tasks/results/star_rsem/kidney_14.5_RNA_1.genes.results",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
annotation <- as.data.frame(read.table("/vol/COMPEPIWS/pipelines/references/mm10_genome_genes.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed6_file <- merge(gene_result, annotation, by.x = "transcript_id.s.", by.y = "V4")
bed6_file <- data.frame(bed6_file$V1, bed6_file$V2, bed6_file$V3, bed6_file$gene_id, bed6_file$TPM, bed6_file$V6)
colnames(bed6_file) <- c("chr", "start", "end", "gene_ID", "TPM", "strand")

write.table(bed6_file, "kidney_14.5_rep1.bed", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

```

``` bash
# b)
# kidney_14.5_H3K27ac_R1
$ computeMatrix scale-regions -S /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K27ac_R1.bigWig  
                                /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K27me3_R1.bigWig 
                                /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K36me3_R1.bigWig 
                                /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me1_R1.bigWig
                                /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me3_R1.bigWig 
                                /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K9ac_R1.bigWig 
                                /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K9me3_R1.bigWig          
                              -R /vol/COMPEPIWS/groups/rnaseq1/tasks/Siwen_work/kidney_14.5_rep1.bed 
                              --beforeRegionStartLength 3000 
                              --regionBodyLength 5000 
                              --afterRegionStartLength 3000 
                              --skipZeros 
                              -o matrix.mat.gz   
                            
$ plotHeatmap -m matrix_all.mat.gz -out Heatmap_seven_histone.png 
```

Here we plot a overall heatmap, also a single example for one histone mark.

![](https://i.imgur.com/auYprdW.jpg)

![](https://i.imgur.com/7C8SFo7.png)


c) Marks H3K27ac, H3K4me3 and H3K9ac are enriched ar the TSS, and mark H3K36me3 is entiched in the gene body. 

d) 

``` bash
# e)
$ plotHeatmap -m matrix_all.mat.gz 
              -out kmeans_all.png 
              --colorMap RdBu 
              --whatToShow 'heatmap and colorbar' 
              --kmeans 6
              --outFileSortedRegions Heatmap1sortedRegions.bed

```

![](https://i.imgur.com/XnkIuD8.jpg)

![](https://i.imgur.com/DmpYob4.png)

``` R
# d)
kmeans_clusters <- as.data.frame(read.table("/vol/COMPEPIWS/groups/rnaseq1/tasks/Siwen_work/results/IntegratedAnslysis-6/kmeans_cluster.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
df <- data.frame(kmeans_clusters$V5, kmeans_clusters$V13)
colnames(df) <- c("TPM", "cluster")
png("boxplot_TPM.PNG", width=16,height=10,units="in",res=600)
par(mar=c(7.1, 7.1, 7.1, 3.1))
colnames(df) <- c("TPM","cluster")
ggplot(df, aes(x=cluster, y=TPM)) + geom_boxplot()
dev.off()
 
```

![](https://i.imgur.com/ClNLg42.png)


### 7 Chromatin accessibility and gene expression

``` R
# a.
peaks_df = read.table("/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq_ref/differential_peaks.tsv", header=TRUE)
degs_complete = read.table("/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq1/DEGs/DEGs_complete.bed", header=TRUE, sep = "\t")

peaks_df_common = match(degs_complete$nearestGene, peaks_df$gene_ID, nomatch="") 
peaks_df_common = unique(peaks_df_common[!is.na(peaks_df_common)])
peaks_df_common = peaks_df[peaks_df_common, ]
peaks_df_common = peaks_df_common[order(peaks_df_common$gene_name), ]

degs_complete_common = match(peaks_df$gene_name, degs_complete$gene_names, nomatch="") 
degs_complete_common = unique(degs_complete_common[!is.na(degs_complete_common)])
degs_complete_common = degs_complete[degs_complete_common, ]
degs_complete_common = degs_complete_common[order(degs_complete_common$gene_names), ]

combined_table = merge(degs_complete_common, peaks_df_common, by.x="gene_names", by.y="gene_name")
deg_diff_peak_combine = data.frame(combined_table$gene_names, combined_table$distance, combined_table$meanCPMLiver, combined_table$meanCPMKidney, combined_table$meanLog10FpkmGrp2_liver, combined_table$meanLog10FpkmGrp1_kidney, combined_table$logFC, combined_table$log2FoldChange)
colnames(combined_table) = c("gene_name", "distance", "meanCPMLiver", "meanCPMKidney", "meanLog10FpkmGrp2_liver", "meanLog10FpkmGrp1_kidney", "logFC", "log2FoldChange")

```

### 9 Enhancers and gene expression 

``` R
# a)
library(data.table)

enh_kid <- read.table("/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/segmentation/kidney_15_segments_dense.bed", head=FALSE, sep="\t", fill=TRUE)
enh_kid <- as.data.frame(enh_kid)
enh_kid <- enh_kid[-1,]
enh_kid <- enh_kid[enh_kid$V4 %like% "Enh", ]
enh_kid <- data.frame(enh_kid$V1, enh_kid$V2, enh_kid$V3)
colnames(enh_kid) <- c("chrom", "chromStart", "chromEnd")
write.table(enh_kid, "enh_kid.bed", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


enh_liver <- read.table("/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/segmentation/liver_15_segments_dense.bed", head=FALSE, sep="\t", fill=TRUE)
enh_liver <- as.data.frame(enh_liver)
enh_liver <- enh_liver[-1,]
enh_liver <- enh_liver[enh_liver$V4 %like% "Enh", ]
enh_liver <- data.frame(enh_liver$V1, enh_liver$V2, enh_liver$V3)
colnames(enh_liver) <- c("chrom", "chromStart", "chromEnd")
write.table(enh_liver, "enh_liver.bed", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
```

``` bash
# b)
$ bedtools merge -i enh_kid.bed > enh_kid_merge.bed
$ bedtools merge -i enh_liver.bed > enh_liver_merge.bed
```

``` R
# d)




```






