# <!------------- New analysis from the beginning ----------!> 

# Codes that are only needed once is folded in the sections to save the space and make it easier for double check

# load pre-treated data or raw counts, change the names of variables -------------------------------------


setwd('D:/foldername/')
print('Set a debug point to stop accidently ctrl+shift+enter')
# load the mouse Cell in vitro 2D data from Lianne
load(file = 'Cell_Dedif/20210817_Lianne_Cell_compound_2D.RData')
dds_2D_mCell <- RNA_Seq_Cell2  # wrong name of dds, change it
rm(RNA_Seq_Cell2)
rld_2D_mCell <- rld
rm(rld)
rm(Annotation)  # convert annotation by myself (avoid old version of annotation)


# load the mouse Cell in vivo data 
load(file = 'Cell_Dedif/foldername_Seq_acute_chronic_Treat_Cell.RData')
dds_invivo_mCell <- RNA_Seq_Cell2  # wrong name of dds, change it
rm(RNA_Seq_Cell2)
rld_invivo_mCell <- rld
rm(rld)


# load the mouse co-culture 3Ds data (D0 - D10)
load(file = 'Mouse_3D_timeline/ImpulseDE2_script/annot_all_matched.RData')
load(file = 'Mouse_3D_timeline/ImpulseDE2_script/counts_3Ds.RData')
head(counts_3Ds)



# Create in vitro 2D culture counts, dds, rld and annot for Wald test -----------------------------------------------
library(DESeq2)
dds_2D_mCell
counts_2D_mCell <- counts(dds_2D_mCell)  # get the raw counts
head(counts_2D_mCell)

# remove very low counts and replace the rownames of counts into gene symbol
# there are more than 10 samples in the count matrix, so if rowsums less than 10, it means only 1 or 2 count in 1 condition
counts_2D_mCell_symbol <- counts_2D_mCell
counts_2D_mCell_symbol <- counts_2D_mCell_symbol[rowSums(counts_2D_mCell_symbol)>10, ]  
nrow(counts_2D_mCell_symbol)
nrow(counts_2D_mCell)

# start convert ensembl id to gene symbols
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

genelist_counts_2D_mCell_symbol <- rownames(counts_2D_mCell_symbol)
head(genelist_counts_2D_mCell_symbol)
length(genelist_counts_2D_mCell_symbol)

annot_2D <- getBM(attributes = c(
  'mgi_symbol', 'ensembl_gene_id'), 
  filters = 'ensembl_gene_id', 
  values=genelist_counts_2D_mCell_symbol,
  mart = ensembl
)

head(annot_2D)

# Select the rownames (geneid) did not be matched by biomaRt
mannually_annot_2D <- subset(rownames(counts_2D_mCell_symbol), !(rownames(counts_2D_mCell_symbol) %in% annot_2D$ensembl_gene_id))
head(mannually_annot_2D)
length(mannually_annot_2D)
mannually_annot_2D <- data.frame(t(mannually_annot_2D))
mannually_annot_2D <- t(mannually_annot_2D)
head(mannually_annot_2D)
# View(mannually_annot)
nrow(mannually_annot_2D)
setwd('D:/foldername/Cell_Dedif/20230906/Temporary files from R/')
write.csv(mannually_annot_2D, 'unmatched annot 2D.csv', row.names = FALSE)
# after https://www.biotools.fr/mouse/ensembl_symbol_converter convert the unmatched annot, read the csv again
mannually_annot_2D_finished <- read.csv('unmatched annot 2D.csv', header = TRUE, sep = ';')
colnames(mannually_annot_2D_finished) <- c('ensembl_gene_id', 'mgi_symbol')
head(mannually_annot_2D_finished)
nrow(mannually_annot_2D_finished)
tail(mannually_annot_2D_finished)

annot_2D_all_matched <- rbind(annot_2D, mannually_annot_2D_finished)
head(annot_2D_all_matched)
tail(annot_2D_all_matched)
# check how many rows in mgi_symbol are empty, 
# and remove them or replace them by https://hostdb.org/hostdb/app/search?q=ENSMUSG00000094915
sum(annot_2D_all_matched$mgi_symbol == '')
View(annot_2D_all_matched[annot_2D_all_matched$mgi_symbol == '', ])
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000062783', ]$mgi_symbol <- 'Csprs'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000079222', ]$mgi_symbol <- 'AC132444.1'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000079800', ]$mgi_symbol <- 'AC125149.3'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000079808', ]$mgi_symbol <- 'AC168977.1'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000083658', ]$mgi_symbol <- 'Gm15798'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000094054', ]$mgi_symbol <- 'AC133103.2'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000094514', ]$mgi_symbol <- 'AC133103.3'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000094915', ]$mgi_symbol <- 'AC168977.2'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000095041', ]$mgi_symbol <- 'AC149090.1'
annot_2D_all_matched[annot_2D_all_matched$ensembl_gene_id == 'ENSMUSG00000115902', ]$mgi_symbol <- 'AC113595.1'

annot_2D_all_matched <- annot_2D_all_matched[!(annot_2D_all_matched$mgi_symbol == ''), ]
head(annot_2D_all_matched)
sum(annot_2D_all_matched$mgi_symbol == '')
nrow(annot_2D_all_matched)

# check completely duplicated rows
sum(duplicated(annot_2D_all_matched)) 

# sort annot as the same order as in the counts (ensembl gene id order)
annot_2D_all_matched_sorted <- annot_2D_all_matched[order(annot_2D_all_matched$ensembl_gene_id), ]
head(annot_2D_all_matched_sorted)
sum(!(annot_2D_all_matched_sorted$ensembl_gene_id == rownames(counts_2D_mCell_symbol))) # check if the order matched

# Check duplicated mgi_symbol :(
sum(duplicated(annot_2D_all_matched_sorted$mgi_symbol))
duplicated_2D_mgi_symbol <- annot_2D_all_matched_sorted[
  annot_2D_all_matched_sorted$mgi_symbol %in% 
    unique(annot_2D_all_matched_sorted$mgi_symbol[duplicated(annot_2D_all_matched_sorted$mgi_symbol)]), ]
head(duplicated_2D_mgi_symbol)
nrow(duplicated_2D_mgi_symbol)

# discard the duplicated rownames that are not the main gene transcription
nrow(annot_2D_all_matched_sorted) - nrow(annot_2D_all_matched_sorted[!(annot_2D_all_matched_sorted$ensembl_gene_id %in% 
                                                                   c('ENSMUSG00000115018', 'ENSMUSG00000107877', 
                                                                     'ENSMUSG00000079737', 'ENSMUSG00000091071', 'ENSMUSG00000083012', 
                                                                     'ENSMUSG00000087014', 'ENSMUSG00000097823', 
                                                                     'ENSMUSG00000056089', 'ENSMUSG00000105061', 
                                                                     'ENSMUSG00000109685', 'ENSMUSG00000089945', 
                                                                     'ENSMUSG00000090053', 'ENSMUSG00000026064', 
                                                                     'ENSMUSG00000115420', 'ENSMUSG00000116048', 
                                                                     'ENSMUSG00000116184')), ])
annot_2D_all_matched_sorted <- annot_2D_all_matched_sorted[!(annot_2D_all_matched_sorted$ensembl_gene_id %in% 
                                                         c('ENSMUSG00000115018', 'ENSMUSG00000107877', 
                                                           'ENSMUSG00000079737', 'ENSMUSG00000091071', 'ENSMUSG00000083012', 
                                                           'ENSMUSG00000087014', 'ENSMUSG00000097823', 
                                                           'ENSMUSG00000056089', 'ENSMUSG00000105061', 
                                                           'ENSMUSG00000109685', 'ENSMUSG00000089945', 
                                                           'ENSMUSG00000090053', 'ENSMUSG00000026064', 
                                                           'ENSMUSG00000115420', 'ENSMUSG00000116048', 
                                                           'ENSMUSG00000116184')), ]
nrow(annot_2D_all_matched_sorted)
head(annot_2D_all_matched_sorted)

counts_2D_mCell_symbol_selected <- counts_2D_mCell_symbol[!(rownames(counts_2D_mCell_symbol) %in% 
                                                     c('ENSMUSG00000115018', 'ENSMUSG00000107877', 
                                                       'ENSMUSG00000079737', 'ENSMUSG00000091071', 'ENSMUSG00000083012', 
                                                       'ENSMUSG00000087014', 'ENSMUSG00000097823', 
                                                       'ENSMUSG00000056089', 'ENSMUSG00000105061', 
                                                       'ENSMUSG00000109685', 'ENSMUSG00000089945', 
                                                       'ENSMUSG00000090053', 'ENSMUSG00000026064', 
                                                       'ENSMUSG00000115420', 'ENSMUSG00000116048', 
                                                       'ENSMUSG00000116184')), ]
nrow(counts_2D_mCell_symbol_selected)
head(counts_2D_mCell_symbol_selected)

rownames(counts_2D_mCell_symbol_selected) <- annot_2D_all_matched_sorted$mgi_symbol
head(counts_2D_mCell_symbol_selected)

# The colnames order is terrible, so change them. Never seen such an anti-humanity order. 
counts_2D_mCell_symbol_selected <- counts_2D_mCell_symbol_selected[,c('0_Solv', '0_Solv.1', '0_Solv.2', '0_Solv.3', 
                                                   '1_Solv', '1_Solv.1', '1_Solv.2', '1_Solv.3',
                                                   '2_Solv', '2_Solv.1', '2_Solv.2', '2_Solv.3',
                                                   '4_Solv', '4_Solv.1', '4_Solv.2', '4_Solv.3',
                                                   '1_compound', '1_compound.1', '1_compound.2', '1_compound.3', 
                                                   '2_compound', '2_compound.1', '2_compound.2', '2_compound.3', 
                                                   '4_compound', '4_compound.1', '4_compound.2', '4_compound.3'
                                                   )]
colnames(counts_2D_mCell_symbol_selected) <- c('D0_Solv', 'D0_Solv.1', 'D0_Solv.2', 'D0_Solv.3', 
                                      'D1_Solv', 'D1_Solv.1', 'D1_Solv.2', 'D1_Solv.3',
                                      'D2_Solv', 'D2_Solv.1', 'D2_Solv.2', 'D2_Solv.3',
                                      'D4_Solv', 'D4_Solv.1', 'D4_Solv.2', 'D4_Solv.3',
                                      'D1_compound', 'D1_compound.1', 'D1_compound.2', 'D1_compound.3', 
                                      'D2_compound', 'D2_compound.1', 'D2_compound.2', 'D2_compound.3', 
                                      'D4_compound', 'D4_compound.1', 'D4_compound.2', 'D4_compound.3')
nrow(counts_2D_mCell_symbol_selected)


# Start DESeq2
samplelist_2D_mCell <- data.frame(Timepoint = factor(c(rep(c('0'), 4), rep(c('1'),4), 
                                                       rep(c('2'), 4), rep(c('4'), 4),
                                                       rep(c('1'),4), rep(c('2'), 4), rep(c('4'), 4)
                                                       )
                                                     ),
                                  Treatment = factor(c(rep(c('Solv'), 16), rep(c('compound'), 12)
                                                      )
                                                    )
                                  )
rownames(samplelist_2D_mCell) <- colnames(counts_2D_mCell_symbol_selected)
head(samplelist_2D_mCell)
# dds_2D_mCell <- DESeqDataSetFromMatrix(countData = as.matrix(counts_2D_mCell_symbol_selected), 
#                                        colData = samplelist_2D_mCell, 
#                                        design = ~ Treatment + Timepoint + Treatment:Timepoint)


# here a error given: 
  # Error in checkFullRank(modelMatrix) : 
  #   the model matrix is not full rank, so the model cannot be fit as specified.
  # One or more variables or interaction terms in the design formula are linear
  # combinations of the others and must be removed.
  # 
  # Please read the vignette section 'Model matrix not full rank':
  #   
  #   vignette('DESeq2')
  # In addition: Warning message:
  #   In DESeqDataSet(se, design = design, ignoreRank) :
  #   1 duplicate rownames were renamed by adding numbers


# caused by the Timepoint = 0 with treatment ctrl/solv, so it is better to first see the timepoint change at the D0 and D4 in 2 groups

# <!----- Here is the dds of timepoint D0 and D4--------------!>
dds_2D_mCell_timepoint <- DESeqDataSetFromMatrix(countData = as.matrix(counts_2D_mCell_symbol_selected), 
                                                 colData = samplelist_2D_mCell, 
                                                 design = ~ Timepoint)

dds_2D_mCell_timepoint <- DESeq(dds_2D_mCell_timepoint)

# deseq2 choose a reference level for factors based on alphabetical order, so have to tell deseq2 which one should be compared
dds_2D_mCell_timepoint$Timepoint <- relevel(dds_2D_mCell_timepoint$Timepoint, ref = '0')  # is it necessary to claim the ref before subset?
rld_2D_mCell_timepoint <- rlog(dds_2D_mCell_timepoint)
normcounts_2D_mCell_timepoint <- counts(dds_2D_mCell_timepoint, normalized = TRUE)

# subset the solv/compound treatment and drop the levels that don't have the samples, then specific the contrast
dds_2D_mCell_timepoint_solv <- dds_2D_mCell_timepoint[, dds_2D_mCell_timepoint$Treatment == 'Solv']
dds_2D_mCell_timepoint_solv$Treatment <- droplevels(dds_2D_mCell_timepoint_solv$Treatment)
dds_2D_mCell_timepoint_solv_D4 <- dds_2D_mCell_timepoint_solv
dds_2D_mCell_timepoint_solv_D4$Timepoint <- factor(dds_2D_mCell_timepoint_solv_D4$Timepoint, levels = c('0', '4'))
res_2D_mCell_timepoint_solv_D4 <- results(dds_2D_mCell_timepoint_solv_D4, alpha = 0.05, contrast = c('Timepoint', '4','0'))   
# alpha for padjust cutoff (FDR), contrast again
summary(res_2D_mCell_timepoint_solv_D4)
rld_2D_mCell_timepoint_solv_D4 <- rlog(dds_2D_mCell_timepoint_solv_D4) 
normcounts_2D_mCell_timepoint_solv_D4 <- counts(dds_2D_mCell_timepoint_solv_D4, normalized= TRUE)
# tired of so many different groups. The names of variables are as long as the toilet paper


dds_2D_mCell_timepoint_compound <- dds_2D_mCell_timepoint[, dds_2D_mCell_timepoint$Treatment == 'compound' | 
                                                        dds_2D_mCell_timepoint$Timepoint == '0']
dds_2D_mCell_timepoint_compound$Treatment <- droplevels(dds_2D_mCell_timepoint_compound$Treatment)
dds_2D_mCell_timepoint_compound_D4 <- dds_2D_mCell_timepoint_compound
dds_2D_mCell_timepoint_compound_D4$Timepoint <- factor(dds_2D_mCell_timepoint_compound_D4$Timepoint, levels = c('0', '4'))
res_2D_mCell_timepoint_compound_D4 <- results(dds_2D_mCell_timepoint_compound_D4, alpha = 0.05, contrast = c('Timepoint', '4','0')) # alpha for padjust cutoff (FDR)
summary(res_2D_mCell_timepoint_compound_D4)
rld_2D_mCell_timepoint_compound_D4 <- rlog(dds_2D_mCell_timepoint_compound_D4)
normcounts_2D_mCell_timepoint_compound_D4 <- counts(dds_2D_mCell_timepoint_compound_D4, normalized = TRUE)
# DO NOT try to change these variables when haven't waken up. Can never figure out them without a clear mind



# <!----- Here is the dds of timepoint D0 and D1--------------!>

# subset the solv/compound treatment and drop the levels that don't have the samples, then specific the contrast
# dds_2D_mCell_timepoint_solv <- dds_2D_mCell_timepoint[, dds_2D_mCell_timepoint$Treatment == 'Solv']
# dds_2D_mCell_timepoint_solv$Treatment <- droplevels(dds_2D_mCell_timepoint_solv$Treatment)
dds_2D_mCell_timepoint_solv_D1 <- dds_2D_mCell_timepoint_solv
dds_2D_mCell_timepoint_solv_D1$Timepoint <- factor(dds_2D_mCell_timepoint_solv_D1$Timepoint, levels = c('0', '1'))
res_2D_mCell_timepoint_solv_D1 <- results(dds_2D_mCell_timepoint_solv_D1, alpha = 0.05, contrast = c('Timepoint', '1','0'))   # alpha for padjust cutoff (FDR)
summary(res_2D_mCell_timepoint_solv_D1)
rld_2D_mCell_timepoint_solv_D1 <- rlog(dds_2D_mCell_timepoint_solv_D1) 
normcounts_2D_mCell_timepoint_solv_D1 <- counts(dds_2D_mCell_timepoint_solv_D1, normalized= TRUE)
# tired of so many different groups. The names of variables are as long as the toilet paper


# dds_2D_mCell_timepoint_compound <- dds_2D_mCell_timepoint[, dds_2D_mCell_timepoint$Treatment == 'compound']
# dds_2D_mCell_timepoint_compound$Treatment <- droplevels(dds_2D_mCell_timepoint_compound$Treatment)
dds_2D_mCell_timepoint_compound_D1 <- dds_2D_mCell_timepoint_compound
dds_2D_mCell_timepoint_compound_D1$Timepoint <- factor(dds_2D_mCell_timepoint_compound_D1$Timepoint, levels = c('0', '1'))
res_2D_mCell_timepoint_compound_D1 <- results(dds_2D_mCell_timepoint_compound_D1, alpha = 0.05, contrast = c('Timepoint', '1','0')) # alpha for padjust cutoff (FDR)
summary(res_2D_mCell_timepoint_compound_D1)
rld_2D_mCell_timepoint_compound_D1 <- rlog(dds_2D_mCell_timepoint_compound_D1)
normcounts_2D_mCell_timepoint_compound_D1 <- counts(dds_2D_mCell_timepoint_compound_D1, normalized = TRUE)
# DO NOT try to change these variables when haven't waken up. Can never figure out them without a clear mind

# <!----- Here is the dds of timepoint D0 and D2--------------!>

# subset the solv/compound treatment and drop the levels that don't have the samples, then specific the contrast
# dds_2D_mCell_timepoint_solv <- dds_2D_mCell_timepoint[, dds_2D_mCell_timepoint$Treatment == 'Solv']
# dds_2D_mCell_timepoint_solv$Treatment <- droplevels(dds_2D_mCell_timepoint_solv$Treatment)
dds_2D_mCell_timepoint_solv_D2 <- dds_2D_mCell_timepoint_solv
dds_2D_mCell_timepoint_solv_D2$Timepoint <- factor(dds_2D_mCell_timepoint_solv_D2$Timepoint, levels = c('0', '2'))
res_2D_mCell_timepoint_solv_D2 <- results(dds_2D_mCell_timepoint_solv_D2, alpha = 0.05, contrast = c('Timepoint', '2','0'))   # alpha for padjust cutoff (FDR)
summary(res_2D_mCell_timepoint_solv_D2)
rld_2D_mCell_timepoint_solv_D2 <- rlog(dds_2D_mCell_timepoint_solv_D2) 
normcounts_2D_mCell_timepoint_solv_D2 <- counts(dds_2D_mCell_timepoint_solv_D2, normalized= TRUE)
# tired of so many different groups. The names of variables are as long as the toilet paper


# dds_2D_mCell_timepoint_compound <- dds_2D_mCell_timepoint[, dds_2D_mCell_timepoint$Treatment == 'compound']
# dds_2D_mCell_timepoint_compound$Treatment <- droplevels(dds_2D_mCell_timepoint_compound$Treatment)
dds_2D_mCell_timepoint_compound_D2 <- dds_2D_mCell_timepoint_compound
dds_2D_mCell_timepoint_compound_D2$Timepoint <- factor(dds_2D_mCell_timepoint_compound_D2$Timepoint, levels = c('0', '2'))
res_2D_mCell_timepoint_compound_D2 <- results(dds_2D_mCell_timepoint_compound_D2, alpha = 0.05, contrast = c('Timepoint', '2','0')) # alpha for padjust cutoff (FDR)
summary(res_2D_mCell_timepoint_compound_D2)
rld_2D_mCell_timepoint_compound_D2 <- rlog(dds_2D_mCell_timepoint_compound_D2)
normcounts_2D_mCell_timepoint_compound_D2 <- counts(dds_2D_mCell_timepoint_compound_D2, normalized = TRUE)
# DO NOT try to change these variables when haven't waken up. Can never figure out them without a clear mind

# Leave the time course/transient changes after these analysis 


# Create in vivo acute injury counts, dds, rld and annot for Wald test -----------------------------------------------
library(DESeq2)
dds_invivo_mCell
counts_invivo_mCell <- counts(dds_invivo_mCell)  # get the raw counts
head(counts_invivo_mCell)

# remove very low counts and replace the rownames of counts into gene symbol
# there are more than 10 samples in the count matrix, so if rowsums less than 10, it means only 1 or 2 count in 1 condition
counts_invivo_mCell_symbol <- counts_invivo_mCell
counts_invivo_mCell_symbol <- counts_invivo_mCell_symbol[rowSums(counts_invivo_mCell_symbol)>10, ]  
nrow(counts_invivo_mCell_symbol)
nrow(counts_invivo_mCell)

# start convert ensembl id to gene symbols
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

genelist_counts_invivo_mCell_symbol <- rownames(counts_invivo_mCell_symbol)
head(genelist_counts_invivo_mCell_symbol)
length(genelist_counts_invivo_mCell_symbol)

annot_invivo <- getBM(attributes = c(
  'mgi_symbol', 'ensembl_gene_id'), 
  filters = 'ensembl_gene_id', 
  values=genelist_counts_invivo_mCell_symbol,
  mart = ensembl
)

head(annot_invivo)

# Select the rownames (geneid) did not be matched by biomaRt
mannually_annot_invivo <- subset(rownames(counts_invivo_mCell_symbol), !(rownames(counts_invivo_mCell_symbol) %in% annot_invivo$ensembl_gene_id))
head(mannually_annot_invivo)
length(mannually_annot_invivo)
mannually_annot_invivo <- data.frame(t(mannually_annot_invivo))
mannually_annot_invivo <- t(mannually_annot_invivo)
head(mannually_annot_invivo)
# View(mannually_annot)
nrow(mannually_annot_invivo)
setwd('D:/foldername/Cell_Dedif/20230906/Temporary files from R/')
write.csv(mannually_annot_invivo, 'unmatched annot invivo.csv', row.names = FALSE)
# after https://www.biotools.fr/mouse/ensembl_symbol_converter convert the unmatched annot, read the csv again
mannually_annot_invivo_finished <- read.csv('unmatched annot invivo.csv', header = TRUE, sep = ';')
colnames(mannually_annot_invivo_finished) <- c('ensembl_gene_id', 'mgi_symbol')
head(mannually_annot_invivo_finished)
nrow(mannually_annot_invivo_finished)
tail(mannually_annot_invivo_finished)

annot_invivo_all_matched <- rbind(annot_invivo, mannually_annot_invivo_finished)
head(annot_invivo_all_matched)
tail(annot_invivo_all_matched)
# check how many rows in mgi_symbol are empty, 
# and remove them or replace them by https://hostdb.org/hostdb/app/search?q=ENSMUSG00000094915
sum(annot_invivo_all_matched$mgi_symbol == '')
View(annot_invivo_all_matched[annot_invivo_all_matched$mgi_symbol == '', ])
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000079800', ]$mgi_symbol <- 'AC125149.3'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000079808', ]$mgi_symbol <- 'AC168977.1'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095742', ]$mgi_symbol <- 'CAAA01147332.1'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000079794', ]$mgi_symbol <- 'AC125149.2'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000094799', ]$mgi_symbol <- 'AC125149.4'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000094874', ]$mgi_symbol <- 'AC132444.3'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095041', ]$mgi_symbol <- 'AC149090.1'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095250', ]$mgi_symbol <- 'AC133103.4'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000108728', ]$mgi_symbol <- 'OR5BS1P'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000079192', ]$mgi_symbol <- 'AC125149.1'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000094054', ]$mgi_symbol <- 'AC133103.2'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095092', ]$mgi_symbol <- 'AC125149.5'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000096100', ]$mgi_symbol <- 'AC133103.7'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000079222', ]$mgi_symbol <- 'AC132444.1'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000083658', ]$mgi_symbol <- 'Gm15798'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000094514', ]$mgi_symbol <- 'AC133103.3'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095500', ]$mgi_symbol <- 'AC132444.5'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095672', ]$mgi_symbol <- 'AC133103.5'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095787', ]$mgi_symbol <- 'AC133103.6'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000062783', ]$mgi_symbol <- 'Csprs'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000094728', ]$mgi_symbol <- 'AC132444.2'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000094915', ]$mgi_symbol <- 'AC168977.2'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000096808', ]$mgi_symbol <- 'AC132444.6'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000079190', ]$mgi_symbol <- 'AC133103.1'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000095450', ]$mgi_symbol <- 'AC132444.4'
annot_invivo_all_matched[annot_invivo_all_matched$ensembl_gene_id == 'ENSMUSG00000115902', ]$mgi_symbol <- 'AC113595.1'


annot_invivo_all_matched <- annot_invivo_all_matched[!(annot_invivo_all_matched$mgi_symbol == ''), ]
head(annot_invivo_all_matched)
sum(annot_invivo_all_matched$mgi_symbol == '')
nrow(annot_invivo_all_matched)

# check completely duplicated rows
sum(duplicated(annot_invivo_all_matched)) 

# sort annot as the same order as in the counts (ensembl gene id order)
annot_invivo_all_matched_sorted <- annot_invivo_all_matched[order(annot_invivo_all_matched$ensembl_gene_id), ]
head(annot_invivo_all_matched_sorted)
sum(!(annot_invivo_all_matched_sorted$ensembl_gene_id == rownames(counts_invivo_mCell_symbol))) # check if the order matched

# Check duplicated mgi_symbol :(
sum(duplicated(annot_invivo_all_matched_sorted$mgi_symbol))
duplicated_invivo_mgi_symbol <- annot_invivo_all_matched_sorted[
  annot_invivo_all_matched_sorted$mgi_symbol %in% 
    unique(annot_invivo_all_matched_sorted$mgi_symbol[duplicated(annot_invivo_all_matched_sorted$mgi_symbol)]), ]
head(duplicated_invivo_mgi_symbol)
nrow(duplicated_invivo_mgi_symbol)
View(duplicated_invivo_mgi_symbol)

# discard the duplicated rownames that are not the main gene transcription
nrow(annot_invivo_all_matched_sorted) - nrow(annot_invivo_all_matched_sorted[!(annot_invivo_all_matched_sorted$ensembl_gene_id %in% 
                                                                         c('ENSMUSG00000065361', 'ENSMUSG00000091071', 
                                                                           'ENSMUSG00000107877', 'ENSMUSG00000114515', 
                                                                           'ENSMUSG00000079737', 'ENSMUSG00000112189',
                                                                           'ENSMUSG00000094661', 'ENSMUSG00000095320', 
                                                                           'ENSMUSG00000096271', 'ENSMUSG00000096873', 
                                                                           'ENSMUSG00000095247', 'ENSMUSG00000116429', 
                                                                           'ENSMUSG00000098326', 'ENSMUSG00000115067', 
                                                                           'ENSMUSG00000087014', 'ENSMUSG00000096975', 
                                                                           'ENSMUSG00000097823', 'ENSMUSG00000088022', 
                                                                           'ENSMUSG00000105884', 'ENSMUSG00000106272', 
                                                                           'ENSMUSG00000100992', 'ENSMUSG00000095395',
                                                                           'ENSMUSG00000105863', 'ENSMUSG00000105448', 
                                                                           'ENSMUSG00000065211', 'ENSMUSG00000105690', 
                                                                           'ENSMUSG00000105120', 'ENSMUSG00000105060', 
                                                                           'ENSMUSG00000113702', 'ENSMUSG00000113781', 
                                                                           'ENSMUSG00000056089', 'ENSMUSG00000105061', 
                                                                           'ENSMUSG00000109685', 'ENSMUSG00000116045', 
                                                                           'ENSMUSG00000095623', 'ENSMUSG00000095456', 
                                                                           'ENSMUSG00000104618', 'ENSMUSG00000105743', 
                                                                           'ENSMUSG00000105428', 'ENSMUSG00000104755', 
                                                                           'ENSMUSG00000092920', 'ENSMUSG00000114429', 
                                                                           'ENSMUSG00000115018', 'ENSMUSG00000115074', 
                                                                           'ENSMUSG00000116207', 'ENSMUSG00000089945', 
                                                                           'ENSMUSG00000090053', 'ENSMUSG00000026064', 
                                                                           'ENSMUSG00000115420', 'ENSMUSG00000098943', 
                                                                           'ENSMUSG00000105115', 'ENSMUSG00000099291', 
                                                                           'ENSMUSG00000104896', 'ENSMUSG00000112285', 
                                                                           'ENSMUSG00000116048', 'ENSMUSG00000097052', 
                                                                           'ENSMUSG00000065735', 'ENSMUSG00000110170', 
                                                                           'ENSMUSG00000109224', 'ENSMUSG00000102049',
                                                                           'ENSMUSG00000116275', 'ENSMUSG00000116184'
                                                                           )), ])
annot_invivo_all_matched_sorted <- annot_invivo_all_matched_sorted[!(annot_invivo_all_matched_sorted$ensembl_gene_id %in% 
                                                                       c('ENSMUSG00000065361', 'ENSMUSG00000091071', 
                                                                         'ENSMUSG00000107877', 'ENSMUSG00000114515', 
                                                                         'ENSMUSG00000079737', 'ENSMUSG00000112189',
                                                                         'ENSMUSG00000094661', 'ENSMUSG00000095320', 
                                                                         'ENSMUSG00000096271', 'ENSMUSG00000096873', 
                                                                         'ENSMUSG00000095247', 'ENSMUSG00000116429', 
                                                                         'ENSMUSG00000098326', 'ENSMUSG00000115067', 
                                                                         'ENSMUSG00000087014', 'ENSMUSG00000096975', 
                                                                         'ENSMUSG00000097823', 'ENSMUSG00000088022', 
                                                                         'ENSMUSG00000105884', 'ENSMUSG00000106272', 
                                                                         'ENSMUSG00000100992', 'ENSMUSG00000095395',
                                                                         'ENSMUSG00000105863', 'ENSMUSG00000105448', 
                                                                         'ENSMUSG00000065211', 'ENSMUSG00000105690', 
                                                                         'ENSMUSG00000105120', 'ENSMUSG00000105060', 
                                                                         'ENSMUSG00000113702', 'ENSMUSG00000113781', 
                                                                         'ENSMUSG00000056089', 'ENSMUSG00000105061', 
                                                                         'ENSMUSG00000109685', 'ENSMUSG00000116045', 
                                                                         'ENSMUSG00000095623', 'ENSMUSG00000095456', 
                                                                         'ENSMUSG00000104618', 'ENSMUSG00000105743', 
                                                                         'ENSMUSG00000105428', 'ENSMUSG00000104755', 
                                                                         'ENSMUSG00000092920', 'ENSMUSG00000114429', 
                                                                         'ENSMUSG00000115018', 'ENSMUSG00000115074', 
                                                                         'ENSMUSG00000116207', 'ENSMUSG00000089945', 
                                                                         'ENSMUSG00000090053', 'ENSMUSG00000026064', 
                                                                         'ENSMUSG00000115420', 'ENSMUSG00000098943', 
                                                                         'ENSMUSG00000105115', 'ENSMUSG00000099291', 
                                                                         'ENSMUSG00000104896', 'ENSMUSG00000112285', 
                                                                         'ENSMUSG00000116048', 'ENSMUSG00000097052', 
                                                                         'ENSMUSG00000065735', 'ENSMUSG00000110170', 
                                                                         'ENSMUSG00000109224', 'ENSMUSG00000102049',
                                                                         'ENSMUSG00000116275', 'ENSMUSG00000116184'
                                                                         )), ]
nrow(annot_invivo_all_matched_sorted)
head(annot_invivo_all_matched_sorted)

counts_invivo_mCell_symbol_selected <- counts_invivo_mCell_symbol[!(rownames(counts_invivo_mCell_symbol) %in% 
                                                                c('ENSMUSG00000065361', 'ENSMUSG00000091071', 
                                                                  'ENSMUSG00000107877', 'ENSMUSG00000114515', 
                                                                  'ENSMUSG00000079737', 'ENSMUSG00000112189',
                                                                  'ENSMUSG00000094661', 'ENSMUSG00000095320', 
                                                                  'ENSMUSG00000096271', 'ENSMUSG00000096873', 
                                                                  'ENSMUSG00000095247', 'ENSMUSG00000116429', 
                                                                  'ENSMUSG00000098326', 'ENSMUSG00000115067', 
                                                                  'ENSMUSG00000087014', 'ENSMUSG00000096975', 
                                                                  'ENSMUSG00000097823', 'ENSMUSG00000088022', 
                                                                  'ENSMUSG00000105884', 'ENSMUSG00000106272', 
                                                                  'ENSMUSG00000100992', 'ENSMUSG00000095395',
                                                                  'ENSMUSG00000105863', 'ENSMUSG00000105448', 
                                                                  'ENSMUSG00000065211', 'ENSMUSG00000105690', 
                                                                  'ENSMUSG00000105120', 'ENSMUSG00000105060', 
                                                                  'ENSMUSG00000113702', 'ENSMUSG00000113781', 
                                                                  'ENSMUSG00000056089', 'ENSMUSG00000105061', 
                                                                  'ENSMUSG00000109685', 'ENSMUSG00000116045', 
                                                                  'ENSMUSG00000095623', 'ENSMUSG00000095456', 
                                                                  'ENSMUSG00000104618', 'ENSMUSG00000105743', 
                                                                  'ENSMUSG00000105428', 'ENSMUSG00000104755', 
                                                                  'ENSMUSG00000092920', 'ENSMUSG00000114429', 
                                                                  'ENSMUSG00000115018', 'ENSMUSG00000115074', 
                                                                  'ENSMUSG00000116207', 'ENSMUSG00000089945', 
                                                                  'ENSMUSG00000090053', 'ENSMUSG00000026064', 
                                                                  'ENSMUSG00000115420', 'ENSMUSG00000098943', 
                                                                  'ENSMUSG00000105115', 'ENSMUSG00000099291', 
                                                                  'ENSMUSG00000104896', 'ENSMUSG00000112285', 
                                                                  'ENSMUSG00000116048', 'ENSMUSG00000097052', 
                                                                  'ENSMUSG00000065735', 'ENSMUSG00000110170', 
                                                                  'ENSMUSG00000109224', 'ENSMUSG00000102049',
                                                                  'ENSMUSG00000116275', 'ENSMUSG00000116184'
                                                              )), ]
nrow(counts_invivo_mCell_symbol_selected)
head(counts_invivo_mCell_symbol_selected)
sum(duplicated(annot_invivo_all_matched_sorted$mgi_symbol))

# sort counts order
counts_invivo_mCell_symbol_selected <- counts_invivo_mCell_symbol_selected[order(rownames(counts_invivo_mCell_symbol_selected)), ]
rownames(counts_invivo_mCell_symbol_selected) <- annot_invivo_all_matched_sorted$mgi_symbol
head(counts_invivo_mCell_symbol_selected)

# The colnames order is terrible, so change them. 
colnames(counts_invivo_mCell_symbol_selected) <- c('Healthy', 'Healthy.1', 'Healthy.2', 'Healthy.3', 
                                                   'D1_Injury', 'D1_Injury.1', 
                                                   'D3_Injury', 'D3_Injury.1', 'D3_Injury.2', 'D3_Injury.3',
                                                   'D7_Injury', 'D7_Injury.1', 'D7_Injury.2', 'D7_Injury.3')
head(counts_invivo_mCell_symbol_selected)


# Start DESeq2
samplelist_invivo_mCell <- data.frame(Timepoint = factor(c(rep(c('0'), 4), rep(c('1'),2), 
                                                         rep(c('3'), 4), rep(c('7'), 4))
                                                         ),
                                      Treatment = factor(c(rep(c('Healthy'), 4), rep(c('Acute injury'), 10)))
                                      )
rownames(samplelist_invivo_mCell) <- colnames(counts_invivo_mCell_symbol_selected)
head(samplelist_invivo_mCell)

# In this case, only time points are the variable. 
# Check 2 different comparison - D1 to D0 to see what happened after acute injury. D7 to D0 to see which kind of genes/pathways did not recover. 

# <!----- Here is the dds--------------!>
dds_invivo_mCell_timepoint <- DESeqDataSetFromMatrix(countData = as.matrix(counts_invivo_mCell_symbol_selected), 
                                                     colData = samplelist_invivo_mCell, 
                                                     design = ~ Timepoint)

dds_invivo_mCell_timepoint <- DESeq(dds_invivo_mCell_timepoint)

# deseq2 choose a reference level for factors based on alphabetical order, so have to tell deseq2 which one should be compared
dds_invivo_mCell_timepoint$Timepoint <- relevel(dds_invivo_mCell_timepoint$Timepoint, ref = '0')  # is it necessary to claim the ref before subset?

# Specific the contrast - D0 and D1
dds_invivo_mCell_timepoint_acute <- dds_invivo_mCell_timepoint
dds_invivo_mCell_timepoint_acute$Timepoint <- factor(dds_invivo_mCell_timepoint_acute$Timepoint, levels = c('0', '1'))
res_invivo_mCell_timepoint_acute <- results(dds_invivo_mCell_timepoint_acute, alpha = 0.05, contrast = c('Timepoint', '1','0'))   # alpha for padjust cutoff (FDR)
summary(res_invivo_mCell_timepoint_acute)
rld_invivo_mCell_timepoint_acute <- rlog(dds_invivo_mCell_timepoint_acute) 
normcounts_invivo_mCell_timepoint_acute <- counts(dds_invivo_mCell_timepoint_acute, normalized= TRUE)

# Specific the contrast - D0 and D3
dds_invivo_mCell_timepoint_D3 <- dds_invivo_mCell_timepoint
dds_invivo_mCell_timepoint_D3$Timepoint <- factor(dds_invivo_mCell_timepoint_D3$Timepoint, levels = c('0', '3'))
res_invivo_mCell_timepoint_D3 <- results(dds_invivo_mCell_timepoint_D3, alpha = 0.05, contrast = c('Timepoint', '3','0')) # alpha for padjust cutoff (FDR)
summary(res_invivo_mCell_timepoint_D3)
rld_invivo_mCell_timepoint_D3 <- rlog(dds_invivo_mCell_timepoint_D3)
normcounts_invivo_mCell_timepoint_D3 <- counts(dds_invivo_mCell_timepoint_D3, normalized = TRUE)


# Specific the contrast - D0 and D7
dds_invivo_mCell_timepoint_recover <- dds_invivo_mCell_timepoint
dds_invivo_mCell_timepoint_recover$Timepoint <- factor(dds_invivo_mCell_timepoint_recover$Timepoint, levels = c('0', '7'))
res_invivo_mCell_timepoint_recover <- results(dds_invivo_mCell_timepoint_recover, alpha = 0.05, contrast = c('Timepoint', '7','0')) # alpha for padjust cutoff (FDR)
summary(res_invivo_mCell_timepoint_recover)
rld_invivo_mCell_timepoint_recover <- rlog(dds_invivo_mCell_timepoint_recover)
normcounts_invivo_mCell_timepoint_recover <- counts(dds_invivo_mCell_timepoint_recover, normalized = TRUE)


# Leave the transient changes after these analysis 


# 3Ds counts, dds, rld and annot, wait until the comparison of mCell finished-------------------------------------------
library(DESeq2)

# add annotation for the counts
counts_3Ds_symbol <- counts_3Ds
nrow(counts_3Ds_symbol)
nrow(annot_all_matched)

annot_all_matched_sorted <- annot_all_matched[order(annot_all_matched$ensembl_gene_id), ]
head(annot_all_matched_sorted)
sum(!(annot_all_matched_sorted$ensembl_gene_id == rownames(counts_3Ds_symbol)))
sum(duplicated(annot_all_matched_sorted$mgi_symbol))
duplicated_mgi_symbol <- annot_all_matched_sorted[
  annot_all_matched_sorted$mgi_symbol %in% 
    unique(annot_all_matched_sorted$mgi_symbol[duplicated(annot_all_matched_sorted$mgi_symbol)]), ]
head(duplicated_mgi_symbol)
nrow(duplicated_mgi_symbol)

nrow(annot_all_matched_sorted) - nrow(annot_all_matched_sorted[!(annot_all_matched_sorted$ensembl_gene_id %in% 
                                                                   c('ENSMUSG00000079737', 'ENSMUSG00000090053', 
                                                                     'ENSMUSG00000109685', 'ENSMUSG00000091071', 'ENSMUSG00000112189', 
                                                                     'ENSMUSG00000113662', 'ENSMUSG00000113702')), ])
annot_all_matched_sorted <- annot_all_matched_sorted[!(annot_all_matched_sorted$ensembl_gene_id %in% 
                                                         c('ENSMUSG00000079737', 'ENSMUSG00000090053', 
                                                           'ENSMUSG00000109685', 'ENSMUSG00000091071', 'ENSMUSG00000112189', 
                                                           'ENSMUSG00000113662', 'ENSMUSG00000113702')), ]
nrow(annot_all_matched_sorted)
head(annot_all_matched_sorted)

counts_3Ds_symbol <- counts_3Ds_symbol[!(rownames(counts_3Ds_symbol) %in% 
                                                       c('ENSMUSG00000079737', 'ENSMUSG00000090053', 
                                                         'ENSMUSG00000109685', 'ENSMUSG00000091071', 'ENSMUSG00000112189', 
                                                         'ENSMUSG00000113662', 'ENSMUSG00000113702')), ]
nrow(annot_all_matched_sorted)

rownames(counts_3Ds_symbol) <- annot_all_matched_sorted$mgi_symbol
head(counts_3Ds_symbol)
counts_3Ds_symbol_select <- counts_3Ds_symbol[rowSums(counts_3Ds_symbol) > 10, ]
annot_all_matched_sorted <- annot_all_matched_sorted[annot_all_matched_sorted$mgi_symbol %in% rownames(counts_3Ds_symbol_select), ]

# last check of nrow of counts_3Ds_symbol and annot_all_matched_sorted
nrow(counts_3Ds_symbol_select)
nrow(annot_all_matched_sorted)

# The colnames change. 
colnames(counts_3Ds_symbol_select) <- c('Day0', 'Day0.1', 'Day0.2', 
                                              'Day4', 'Day4.1', 'Day4.2', 
                                              'Day7', 'Day7.1', 'Day7.2', 
                                              'Day10', 'Day10.1', 'Day10.2',
                                              'Day14', 'Day14.1', 'Day14.2')
head(counts_3Ds_symbol_select)


# start Deseq2 
samplelist_3D <- data.frame(Timepoint = factor(c(rep(c('Day0'), 3), rep(c('Day4'),3), 
                                              rep(c('Day7'), 3), 
                                              rep(c('Day10'), 3), rep(c('Day14'), 3) 
                                            ))
                        )
rownames(samplelist_3D) <- colnames(counts_3Ds_symbol_select)
head(samplelist_3D)
tail(samplelist_3D)


dds_3D <- DESeqDataSetFromMatrix(countData = counts_3Ds_symbol_select, 
                              colData = samplelist_3D, 
                              design = ~Timepoint)

dds_3D
dds_3D <- DESeq(dds_3D)
# rm(annot_all_matched)

# deseq2 choose a reference level for factors based on alphabetical order, so have to tell deseq2 which one should be compared
dds_3D$Timepoint <- relevel(dds_3D$Timepoint, ref = 'Day0')

res_3D_D4 <- results(dds_3D, alpha = 0.05, contrast = c('Timepoint', 'Day4','Day0'))  
res_3D_D14 <- results(dds_3D, alpha = 0.05, contrast = c('Timepoint', 'Day14','Day0'))   
# alpha for padjust cutoff (FDR), contrast again
summary(res_3D_D14)
rld_3D <- rlog(dds_3D) 
normcounts_3D <- counts(dds_3D, normalized= TRUE)



# Analysis of in vitro 2D culture, find DE genes and pathways -----------------------------------------

# PCA plot of both treatments and time points 
library(ggplot2)

# create rld for PCA solv
rld_2D_mCell_solv <- rlog(dds_2D_mCell_timepoint_solv) 
pcaData_2D_solv <- plotPCA(rld_2D_mCell_solv, intgroup = c('Timepoint'), returnData = T)
perVar_2D_solv <- round(100 * attr(pcaData_2D_solv, 'percentVar'))

head(pcaData_2D_solv)
# change the timepoint names, make it more professional
pcaData_2D_solv$group <- factor(c(rep('Fresh isolated', 4), rep('Day1', 4), rep('Day2', 4), rep('Day4', 4)))
# adjust the order
pcaData_2D_solv$group <- factor(pcaData_2D_solv$group, levels = c('Fresh isolated', 'Day1', 'Day2', 'Day4'))

ggplot(pcaData_2D_solv, aes(PC1, PC2, color = group)) + 
        geom_point(size = 4) + 
        xlab(paste0("PC1: ", perVar_2D_solv[1], "% variance")) + 
        ylab(paste0("PC2: ", perVar_2D_solv[2], "% variance")) + 
        coord_fixed()+ 
        ggforce::geom_mark_ellipse(aes(fill = group)) + 
        theme(text = element_text(size = 20))  + ggtitle("PCA plot - 2D culture without compound")+
        xlim(-45, 25) +
        ylim(-30, 30) +
        
        theme(
          panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(color = "black") , 
          aspect.ratio = 13/16
        )

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'PCA_2D_mCell_withoutcompound.tiff', 
  dpi = 300
)

# create rld for PCA compound
rld_2D_mCell_compound <- rlog(dds_2D_mCell_timepoint_compound) 
pcaData_2D_compound <- plotPCA(rld_2D_mCell_compound, intgroup = c('Timepoint'), returnData = T)
perVar_2D_compound <- round(100 * attr(pcaData_2D_compound, 'percentVar'))

head(pcaData_2D_compound)
# change the timepoint names, make it more professional
pcaData_2D_compound$group <- factor(c(rep('Fresh isolated', 4), rep('Day1', 4), rep('Day2', 4), rep('Day4', 4)))
# adjust the order
pcaData_2D_compound$group <- factor(pcaData_2D_compound$group, levels = c('Fresh isolated', 'Day1', 'Day2', 'Day4'))

ggplot(pcaData_2D_compound, aes(PC1, PC2, color = group)) + 
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", perVar_2D_compound[1], "% variance")) + 
  ylab(paste0("PC2: ", perVar_2D_compound[2], "% variance")) + 
  coord_fixed()+ 
  ggforce::geom_mark_ellipse(aes(fill = group)) + 
  theme(text = element_text(size = 20))  + ggtitle("PCA plot - 2D culture with compound")+
  xlim(-50, 75) +
  ylim(-30, 50) +
  
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 13/16
  )

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'PCA_2D_mCell_withcompound.tiff', 
  dpi = 300
)

# Create PCA plot for all 2D groups
# rld_2D_mCell_timepoint <- rlog(dds_2D_mCell_timepoint_timepoint) 
pcaData_2D_timepoint <- plotPCA(rld_2D_mCell_timepoint, intgroup = c('Timepoint', 'Treatment'), returnData = T)
perVar_2D_timepoint <- round(100 * attr(pcaData_2D_timepoint, 'percentVar'))

head(pcaData_2D_timepoint)
# change the timepoint names, make it more professional
pcaData_2D_timepoint$Timepoint <- factor(c(rep('Fresh isolated', 4), rep('Day1', 4), rep('Day2', 4), rep('Day4', 4), 
                                           rep('Day1', 4), rep('Day2', 4), rep('Day4', 4)))
# adjust the order
pcaData_2D_timepoint$Timepoint <- factor(pcaData_2D_timepoint$Timepoint, levels = c('Fresh isolated', 'Day1', 'Day2', 'Day4'))
pcaData_2D_timepoint$Timepoint

pcaData_2D_timepoint$Treatment <- factor(c(rep('Fresh isolated', 4), rep('Solv', 12), rep('compound', 12)), 
                                         levels = c('Fresh isolated', 'Solv', 'compound')
                                         )
pcaData_2D_timepoint$Treatment


ggplot(pcaData_2D_timepoint, aes(PC1, PC2, color = Timepoint, shape = Treatment)) + 
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", perVar_2D_timepoint[1], "% variance")) + 
  ylab(paste0("PC2: ", perVar_2D_timepoint[2], "% variance")) + 
  coord_fixed()+ 
  ggforce::geom_mark_ellipse(aes(fill = Timepoint)) + 
  theme(text = element_text(size = 20))  + ggtitle("PCA plot - 2D culture")+
  xlim(-30, 80) +
  ylim(-60, 20) +
  
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 13/16
  )

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'PCA_2D_mCell_alltimepoint.tiff', 
  dpi = 300
)


# Create PCA plot for 2D culture Day1, solv-compound
pcaData_2D_D1 <- plotPCA(rld_2D_mCell_timepoint[, rld_2D_mCell_timepoint$Timepoint %in% c('1')], 
                         intgroup = 'Treatment', returnData = T)
perVar_2D_D1 <- round(100 * attr(pcaData_2D_D1, 'percentVar'))

head(pcaData_2D_D1)

ggplot(pcaData_2D_D1, aes(PC1, PC2, color = Treatment)) + 
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", perVar_2D_timepoint[1], "% variance")) + 
  ylab(paste0("PC2: ", perVar_2D_timepoint[2], "% variance")) + 
  coord_fixed()+ 
  ggforce::geom_mark_ellipse(aes(fill = Treatment)) + 
  theme(text = element_text(size = 20))  + ggtitle("PCA plot - 2D culture - Day1")+
  xlim(-40, 40) +
  ylim(-30, 40) +
  
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 13/16
  )

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'PCA_2D_mCell_D1.tiff', 
  dpi = 300
)


# Create PCA plot for 2D culture Day2, solv-compound
pcaData_2D_D2 <- plotPCA(rld_2D_mCell_timepoint[, rld_2D_mCell_timepoint$Timepoint %in% c('2')], 
                         intgroup = 'Treatment', returnData = T)
perVar_2D_D2 <- round(100 * attr(pcaData_2D_D2, 'percentVar'))

head(pcaData_2D_D2)

ggplot(pcaData_2D_D2, aes(PC1, PC2, color = Treatment)) + 
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", perVar_2D_timepoint[1], "% variance")) + 
  ylab(paste0("PC2: ", perVar_2D_timepoint[2], "% variance")) + 
  coord_fixed()+ 
  ggforce::geom_mark_ellipse(aes(fill = Treatment)) + 
  theme(text = element_text(size = 20))  + ggtitle("PCA plot - 2D culture - Day2")+
  xlim(-35, 80) +
  ylim(-40, 50) +
  
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 13/16
  )

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'PCA_2D_mCell_D2.tiff', 
  dpi = 300
)


# Create PCA plot for 2D culture Day4, solv-compound
pcaData_2D_D4 <- plotPCA(rld_2D_mCell_timepoint[, rld_2D_mCell_timepoint$Timepoint %in% c('4')], 
                         intgroup = 'Treatment', returnData = T)
perVar_2D_D4 <- round(100 * attr(pcaData_2D_D4, 'percentVar'))

head(pcaData_2D_D4)

ggplot(pcaData_2D_D4, aes(PC1, PC2, color = Treatment)) + 
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", perVar_2D_timepoint[1], "% variance")) + 
  ylab(paste0("PC2: ", perVar_2D_timepoint[2], "% variance")) + 
  coord_fixed()+ 
  ggforce::geom_mark_ellipse(aes(fill = Treatment)) + 
  theme(text = element_text(size = 20))  + ggtitle("PCA plot - 2D culture - Day4")+
  xlim(-50, 75) +
  ylim(-60, 30) +
  
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 13/16
  )

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'PCA_2D_mCell_D4.tiff', 
  dpi = 300
)



# we always expect for high variance in the PCA plot. In some cases, it makes sense. 
# However, when the treatment condition added is not strong enough to cause a huge change, the variance is smaller than we think. 
# About the ntop in the plotPCA function: https://www.biostars.org/p/278989/


# find DE genes

## DE genes and GO pathways of solv D1-D0 -----
summary(res_2D_mCell_timepoint_solv_D1)
  # out of 21775 with nonzero total read count
  # adjusted p-value < 0.05
  # LFC > 0 (up)       : 517, 2.4%
  # LFC < 0 (down)     : 529, 2.4%
  # outliers [1]       : 219, 1%
  # low counts [2]     : 8047, 37%
  # (mean count < 9)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results
res_2D_mCell_timepoint_solv_D1 <- as.data.frame(res_2D_mCell_timepoint_solv_D1)
head(res_2D_mCell_timepoint_solv_D1)
sum(is.na(res_2D_mCell_timepoint_solv_D1$padj))
res_2D_mCell_timepoint_solv_D1 <- na.omit(res_2D_mCell_timepoint_solv_D1)

DEgenelist_2D_mCell_timepoint_solv_D1 <- rownames(res_2D_mCell_timepoint_solv_D1[res_2D_mCell_timepoint_solv_D1$padj < 0.05 & 
                                                                                   abs(res_2D_mCell_timepoint_solv_D1$log2FoldChange) >= 1, ])
head(DEgenelist_2D_mCell_timepoint_solv_D1)
length(DEgenelist_2D_mCell_timepoint_solv_D1)

DE_2D_mCell_timepoint_solv_D1 <- res_2D_mCell_timepoint_solv_D1[res_2D_mCell_timepoint_solv_D1$padj < 0.05 & 
                                                                  abs(res_2D_mCell_timepoint_solv_D1$log2FoldChange) >=1, ]
head(DE_2D_mCell_timepoint_solv_D1)
nrow(DE_2D_mCell_timepoint_solv_D1)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
GOpathway_2D_mCell_timepoint_solv_D1 <- enrichGO(gene = DEgenelist_2D_mCell_timepoint_solv_D1, 
                                                 OrgDb = org.Mm.eg.db, 
                                                 ont = 'BP', 
                                                 keyType = 'SYMBOL',
                                                 pAdjustMethod = 'BH', 
                                                 pvalueCutoff = 0.01, 
                                                 qvalueCutoff = 0.05, 
                                                 readable = TRUE
                                                 )
head(GOpathway_2D_mCell_timepoint_solv_D1)
clusterProfiler::dotplot(GOpathway_2D_mCell_timepoint_solv_D1, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - 2D mCell D1 without compound'
                         ) + theme(text = element_text(size = 15), 
                                   legend.text = element_text(size = 15), 
                                   aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_2D_mCell_timepoint_solv_D1 , file = 'GOpathways_2D_mCell_timepoint_solv_D1.xlsx', 
            row.names = TRUE, col.names = TRUE)

                          
# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_2D_mCell_withoutcompound_D1toD0.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_2D_mCell_solv_D1 <- assay(rld_2D_mCell_timepoint_solv_D1)[rownames(assay(rld_2D_mCell_timepoint_solv_D1)) %in% 
                                                                  DEgenelist_2D_mCell_timepoint_solv_D1, 1:8]
head(DE_normcounts_2D_mCell_solv_D1)
nrow(DE_normcounts_2D_mCell_solv_D1)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_2D_mCell_solv_D1 <- as.data.frame(t(scale(t(DE_normcounts_2D_mCell_solv_D1))))
head(scale_DE_normcounts_2D_mCell_solv_D1)


color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_2D_mCell_solv_D1 <- Heatmap(scale_DE_normcounts_2D_mCell_solv_D1,
                                     col = color_ht,
                                     width = unit(10*ncol(scale_DE_normcounts_2D_mCell_solv_D1), 'mm'), 
                                     # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_solv_D1), 'mm'),
                                     name = 'Normalized counts', 
                                     cluster_columns = FALSE,
                                     cluster_rows = TRUE, 
                                     # show_row_names = FALSE,
                                     column_names_side = 'top', 
                                     # rect_gp = gpar(col = 'white', lwd = 1),
                                     # column_split = column_split_overlap,
                                     show_row_names = FALSE,
                                     row_names_max_width = max_text_width(
                                       rownames(scale_DE_normcounts_2D_mCell_solv_D1),
                                       gp = gpar(fontsize = 12))
)
ht_DE_2D_mCell_solv_D1



## DE genes and GO pathways of solv D2-D0 --------
summary(res_2D_mCell_timepoint_solv_D2)
# out of 21775 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 338, 1.6%
# LFC < 0 (down)     : 349, 1.6%
# outliers [1]       : 219, 1%
# low counts [2]     : 7628, 35%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res_2D_mCell_timepoint_solv_D2 <- as.data.frame(res_2D_mCell_timepoint_solv_D2)
head(res_2D_mCell_timepoint_solv_D2)
sum(is.na(res_2D_mCell_timepoint_solv_D2$padj))
res_2D_mCell_timepoint_solv_D2 <- na.omit(res_2D_mCell_timepoint_solv_D2)

DEgenelist_2D_mCell_timepoint_solv_D2 <- rownames(res_2D_mCell_timepoint_solv_D2[res_2D_mCell_timepoint_solv_D2$padj < 0.05 & 
                                                                                   abs(res_2D_mCell_timepoint_solv_D2$log2FoldChange) >= 1, ])
head(DEgenelist_2D_mCell_timepoint_solv_D2)
length(DEgenelist_2D_mCell_timepoint_solv_D2)

DE_2D_mCell_timepoint_solv_D2 <- res_2D_mCell_timepoint_solv_D2[res_2D_mCell_timepoint_solv_D2$padj < 0.05 & 
                                                                  abs(res_2D_mCell_timepoint_solv_D2$log2FoldChange) >=1, ]
head(DE_2D_mCell_timepoint_solv_D2)
nrow(DE_2D_mCell_timepoint_solv_D2)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
GOpathway_2D_mCell_timepoint_solv_D2 <- enrichGO(gene = DEgenelist_2D_mCell_timepoint_solv_D2, 
                                                 OrgDb = org.Mm.eg.db, 
                                                 ont = 'BP', 
                                                 keyType = 'SYMBOL',
                                                 pAdjustMethod = 'BH', 
                                                 pvalueCutoff = 0.01, 
                                                 qvalueCutoff = 0.05, 
                                                 readable = TRUE
                                                ) 
head(GOpathway_2D_mCell_timepoint_solv_D2)
clusterProfiler::dotplot(GOpathway_2D_mCell_timepoint_solv_D2, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - 2D mCell D2 without compound'
                        ) + theme(text = element_text(size = 15), 
                                  legend.text = element_text(size = 15), 
                                  aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_2D_mCell_timepoint_solv_D2 , file = 'GOpathways_2D_mCell_timepoint_solv_D2.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_2D_mCell_withoutcompound_D2toD0.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_2D_mCell_solv_D2 <- assay(rld_2D_mCell_timepoint_solv_D2)[rownames(assay(rld_2D_mCell_timepoint_solv_D2)) %in% 
                                                                          DEgenelist_2D_mCell_timepoint_solv_D2, c(1:4, 9:12)]
head(DE_normcounts_2D_mCell_solv_D2)
nrow(DE_normcounts_2D_mCell_solv_D2)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_2D_mCell_solv_D2 <- as.data.frame(t(scale(t(DE_normcounts_2D_mCell_solv_D2))))
head(scale_DE_normcounts_2D_mCell_solv_D2)


color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_2D_mCell_solv_D2 <- Heatmap(scale_DE_normcounts_2D_mCell_solv_D2,
                                  col = color_ht,
                                  width = unit(10*ncol(scale_DE_normcounts_2D_mCell_solv_D2), 'mm'), 
                                  # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_solv_D2), 'mm'),
                                  name = 'Normalized counts', 
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE, 
                                  # show_row_names = FALSE,
                                  column_names_side = 'top', 
                                  # rect_gp = gpar(col = 'white', lwd = 1),
                                  # column_split = column_split_overlap,
                                  show_row_names = FALSE,
                                  row_names_max_width = max_text_width(
                                    rownames(scale_DE_normcounts_2D_mCell_solv_D2),
                                    gp = gpar(fontsize = 12))
)
ht_DE_2D_mCell_solv_D2


## DE genes and GO pathways of solv D4-D0 ------
summary(res_2D_mCell_timepoint_solv_D4)
  # out of 21775 with nonzero total read count
  # adjusted p-value < 0.05
  # LFC > 0 (up)       : 355, 1.6%
  # LFC < 0 (down)     : 381, 1.7%
  # outliers [1]       : 219, 1%
  # low counts [2]     : 7628, 35%
  # (mean count < 7)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results
res_2D_mCell_timepoint_solv_D4 <- as.data.frame(res_2D_mCell_timepoint_solv_D4)
head(res_2D_mCell_timepoint_solv_D4)
sum(is.na(res_2D_mCell_timepoint_solv_D4$padj))
res_2D_mCell_timepoint_solv_D4 <- na.omit(res_2D_mCell_timepoint_solv_D4)

DEgenelist_2D_mCell_timepoint_solv_D4 <- rownames(res_2D_mCell_timepoint_solv_D4[res_2D_mCell_timepoint_solv_D4$padj < 0.05 & 
                                                                                   abs(res_2D_mCell_timepoint_solv_D4$log2FoldChange) >= 1, ])
head(DEgenelist_2D_mCell_timepoint_solv_D4)
length(DEgenelist_2D_mCell_timepoint_solv_D4)

DE_2D_mCell_timepoint_solv_D4 <- res_2D_mCell_timepoint_solv_D4[res_2D_mCell_timepoint_solv_D4$padj < 0.05 & 
                                                                  abs(res_2D_mCell_timepoint_solv_D4$log2FoldChange) >=1, ]
head(DE_2D_mCell_timepoint_solv_D4)
nrow(DE_2D_mCell_timepoint_solv_D4)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
GOpathway_2D_mCell_timepoint_solv_D4 <- enrichGO(gene = DEgenelist_2D_mCell_timepoint_solv_D4, 
                                                 OrgDb = org.Mm.eg.db, 
                                                 ont = 'BP', 
                                                 keyType = 'SYMBOL',
                                                 pAdjustMethod = 'BH', 
                                                 pvalueCutoff = 0.01, 
                                                 qvalueCutoff = 0.05, 
                                                 readable = TRUE
)
head(GOpathway_2D_mCell_timepoint_solv_D4)
clusterProfiler::dotplot(GOpathway_2D_mCell_timepoint_solv_D4, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - 2D mCell D4 without compound'
                          ) + theme(text = element_text(size = 15), 
                                    legend.text = element_text(size = 15), 
                                    aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_2D_mCell_timepoint_solv_D4 , file = 'GOpathways_2D_mCell_timepoint_solv_D4.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_2D_mCell_withoutcompound_D4toD0.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_2D_mCell_solv_D4 <- assay(rld_2D_mCell_timepoint_solv_D4)[rownames(assay(rld_2D_mCell_timepoint_solv_D4)) %in% 
                                                                          DEgenelist_2D_mCell_timepoint_solv_D4, c(1:4, 13:16)]
head(DE_normcounts_2D_mCell_solv_D4)
nrow(DE_normcounts_2D_mCell_solv_D4)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_2D_mCell_solv_D4 <- as.data.frame(t(scale(t(DE_normcounts_2D_mCell_solv_D4))))
head(scale_DE_normcounts_2D_mCell_solv_D4)


# color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_2D_mCell_solv_D4 <- Heatmap(scale_DE_normcounts_2D_mCell_solv_D4,
                                  col = color_ht,
                                  width = unit(10*ncol(scale_DE_normcounts_2D_mCell_solv_D4), 'mm'), 
                                  # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_solv_D4), 'mm'),
                                  name = 'Normalized counts', 
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE, 
                                  # show_row_names = FALSE,
                                  column_names_side = 'top', 
                                  # rect_gp = gpar(col = 'white', lwd = 1),
                                  # column_split = column_split_overlap,
                                  show_row_names = FALSE,
                                  row_names_max_width = max_text_width(
                                    rownames(scale_DE_normcounts_2D_mCell_solv_D4),
                                    gp = gpar(fontsize = 12))
)
ht_DE_2D_mCell_solv_D4


## DE genes and GO pathways of compound D1-D0 ------
summary(res_2D_mCell_timepoint_compound_D1)
  # out of 21775 with nonzero total read count
  # adjusted p-value < 0.05
  # LFC > 0 (up)       : 518, 2.4%
  # LFC < 0 (down)     : 531, 2.4%
  # outliers [1]       : 218, 1%
  # low counts [2]     : 8048, 37%
  # (mean count < 9)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results
res_2D_mCell_timepoint_compound_D1 <- as.data.frame(res_2D_mCell_timepoint_compound_D1)
head(res_2D_mCell_timepoint_compound_D1)
sum(is.na(res_2D_mCell_timepoint_compound_D1$padj))
res_2D_mCell_timepoint_compound_D1 <- na.omit(res_2D_mCell_timepoint_compound_D1)

DEgenelist_2D_mCell_timepoint_compound_D1 <- rownames(res_2D_mCell_timepoint_compound_D1[res_2D_mCell_timepoint_compound_D1$padj < 0.05 & 
                                                                                   abs(res_2D_mCell_timepoint_compound_D1$log2FoldChange) >= 1, ])
head(DEgenelist_2D_mCell_timepoint_compound_D1)
length(DEgenelist_2D_mCell_timepoint_compound_D1)

DE_2D_mCell_timepoint_compound_D1 <- res_2D_mCell_timepoint_compound_D1[res_2D_mCell_timepoint_compound_D1$padj < 0.05 & 
                                                                  abs(res_2D_mCell_timepoint_compound_D1$log2FoldChange) >=1, ]
head(DE_2D_mCell_timepoint_compound_D1)
nrow(DE_2D_mCell_timepoint_compound_D1)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
GOpathway_2D_mCell_timepoint_compound_D1 <- enrichGO(gene = DEgenelist_2D_mCell_timepoint_compound_D1, 
                                                 OrgDb = org.Mm.eg.db, 
                                                 ont = 'BP', 
                                                 keyType = 'SYMBOL',
                                                 pAdjustMethod = 'BH', 
                                                 pvalueCutoff = 0.01, 
                                                 qvalueCutoff = 0.05, 
                                                 readable = TRUE)
head(GOpathway_2D_mCell_timepoint_compound_D1)
clusterProfiler::dotplot(GOpathway_2D_mCell_timepoint_compound_D1, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - 2D mCell D1 with compound'
                          ) + theme(text = element_text(size = 15), 
                                    legend.text = element_text(size = 15), 
                                    aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_2D_mCell_timepoint_compound_D1 , file = 'GOpathways_2D_mCell_timepoint_compound_D1.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_2D_mCell_withcompound_D1toD0.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_2D_mCell_compound_D1 <- assay(rld_2D_mCell_timepoint_compound_D1)[rownames(assay(rld_2D_mCell_timepoint_compound_D1)) %in% 
                                                                          DEgenelist_2D_mCell_timepoint_compound_D1, 1:8]
head(DE_normcounts_2D_mCell_compound_D1)
nrow(DE_normcounts_2D_mCell_compound_D1)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_2D_mCell_compound_D1 <- as.data.frame(t(scale(t(DE_normcounts_2D_mCell_compound_D1))))
head(scale_DE_normcounts_2D_mCell_compound_D1)


color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_2D_mCell_compound_D1 <- Heatmap(scale_DE_normcounts_2D_mCell_compound_D1,
                                  col = color_ht,
                                  width = unit(10*ncol(scale_DE_normcounts_2D_mCell_compound_D1), 'mm'), 
                                  # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_compound_D1), 'mm'),
                                  name = 'Normalized counts', 
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE, 
                                  # show_row_names = FALSE,
                                  column_names_side = 'top', 
                                  # rect_gp = gpar(col = 'white', lwd = 1),
                                  # column_split = column_split_overlap,
                                  show_row_names = FALSE,
                                  row_names_max_width = max_text_width(
                                    rownames(scale_DE_normcounts_2D_mCell_compound_D1),
                                    gp = gpar(fontsize = 12))
)
ht_DE_2D_mCell_compound_D1

## DE genes and GO pathways of compound D2-D0 ------
summary(res_2D_mCell_timepoint_compound_D2)
# out of 21775 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 337, 1.5%
# LFC < 0 (down)     : 349, 1.6%
# outliers [1]       : 218, 1%
# low counts [2]     : 7629, 35%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
res_2D_mCell_timepoint_compound_D2 <- as.data.frame(res_2D_mCell_timepoint_compound_D2)
head(res_2D_mCell_timepoint_compound_D2)
sum(is.na(res_2D_mCell_timepoint_compound_D2$padj))
res_2D_mCell_timepoint_compound_D2 <- na.omit(res_2D_mCell_timepoint_compound_D2)

DEgenelist_2D_mCell_timepoint_compound_D2 <- rownames(res_2D_mCell_timepoint_compound_D2[res_2D_mCell_timepoint_compound_D2$padj < 0.05 & 
                                                                                   abs(res_2D_mCell_timepoint_compound_D2$log2FoldChange) >= 1, ])
head(DEgenelist_2D_mCell_timepoint_compound_D2)
length(DEgenelist_2D_mCell_timepoint_compound_D2)

DE_2D_mCell_timepoint_compound_D2 <- res_2D_mCell_timepoint_compound_D2[res_2D_mCell_timepoint_compound_D2$padj < 0.05 & 
                                                                  abs(res_2D_mCell_timepoint_compound_D2$log2FoldChange) >=1, ]
head(DE_2D_mCell_timepoint_compound_D2)
nrow(DE_2D_mCell_timepoint_compound_D2)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
GOpathway_2D_mCell_timepoint_compound_D2 <- enrichGO(gene = DEgenelist_2D_mCell_timepoint_compound_D2, 
                                                 OrgDb = org.Mm.eg.db, 
                                                 ont = 'BP', 
                                                 keyType = 'SYMBOL',
                                                 pAdjustMethod = 'BH', 
                                                 pvalueCutoff = 0.01, 
                                                 qvalueCutoff = 0.05, 
                                                 readable = TRUE)
head(GOpathway_2D_mCell_timepoint_compound_D2)
clusterProfiler::dotplot(GOpathway_2D_mCell_timepoint_compound_D2, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - 2D mCell D2 with compound'
                        ) + theme(text = element_text(size = 15), 
                                  legend.text = element_text(size = 15), 
                                  aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_2D_mCell_timepoint_compound_D2 , file = 'GOpathways_2D_mCell_timepoint_compound_D2.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_2D_mCell_withcompound_D2toD0.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_2D_mCell_compound_D2 <- assay(rld_2D_mCell_timepoint_compound_D2)[rownames(assay(rld_2D_mCell_timepoint_compound_D2)) %in% 
                                                                          DEgenelist_2D_mCell_timepoint_compound_D2, c(1:4, 9:12)]
head(DE_normcounts_2D_mCell_compound_D2)
nrow(DE_normcounts_2D_mCell_compound_D2)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_2D_mCell_compound_D2 <- as.data.frame(t(scale(t(DE_normcounts_2D_mCell_compound_D2))))
head(scale_DE_normcounts_2D_mCell_compound_D2)


color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_2D_mCell_compound_D2 <- Heatmap(scale_DE_normcounts_2D_mCell_compound_D2,
                                  col = color_ht,
                                  width = unit(10*ncol(scale_DE_normcounts_2D_mCell_compound_D2), 'mm'), 
                                  # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_compound_D2), 'mm'),
                                  name = 'Normalized counts', 
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE, 
                                  # show_row_names = FALSE,
                                  column_names_side = 'top', 
                                  # rect_gp = gpar(col = 'white', lwd = 1),
                                  # column_split = column_split_overlap,
                                  show_row_names = FALSE,
                                  row_names_max_width = max_text_width(
                                    rownames(scale_DE_normcounts_2D_mCell_compound_D2),
                                    gp = gpar(fontsize = 12))
                                  )
ht_DE_2D_mCell_compound_D2
## DE genes and GO pathways of compound D4-D0 ------
summary(res_2D_mCell_timepoint_compound_D4)
  # out of 21775 with nonzero total read count
  # adjusted p-value < 0.05
  # LFC > 0 (up)       : 354, 1.6%
  # LFC < 0 (down)     : 381, 1.7%
  # outliers [1]       : 218, 1%
  # low counts [2]     : 7629, 35%
  # (mean count < 7)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results
res_2D_mCell_timepoint_compound_D4 <- as.data.frame(res_2D_mCell_timepoint_compound_D4)
head(res_2D_mCell_timepoint_compound_D4)
sum(is.na(res_2D_mCell_timepoint_compound_D4$padj))
res_2D_mCell_timepoint_compound_D4 <- na.omit(res_2D_mCell_timepoint_compound_D4)

DEgenelist_2D_mCell_timepoint_compound_D4 <- rownames(res_2D_mCell_timepoint_compound_D4[res_2D_mCell_timepoint_compound_D4$padj < 0.05 & 
                                                                                   abs(res_2D_mCell_timepoint_compound_D4$log2FoldChange) >= 1, ])
head(DEgenelist_2D_mCell_timepoint_compound_D4)
length(DEgenelist_2D_mCell_timepoint_compound_D4)

DE_2D_mCell_timepoint_compound_D4 <- res_2D_mCell_timepoint_compound_D4[res_2D_mCell_timepoint_compound_D4$padj < 0.05 & 
                                                                  abs(res_2D_mCell_timepoint_compound_D4$log2FoldChange) >=1, ]
head(DE_2D_mCell_timepoint_compound_D4)
nrow(DE_2D_mCell_timepoint_compound_D4)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
GOpathway_2D_mCell_timepoint_compound_D4 <- enrichGO(gene = DEgenelist_2D_mCell_timepoint_compound_D4, 
                                                 OrgDb = org.Mm.eg.db, 
                                                 ont = 'BP', 
                                                 keyType = 'SYMBOL',
                                                 pAdjustMethod = 'BH', 
                                                 pvalueCutoff = 0.01, 
                                                 qvalueCutoff = 0.05, 
                                                 readable = TRUE
)
head(GOpathway_2D_mCell_timepoint_compound_D4)
clusterProfiler::dotplot(GOpathway_2D_mCell_timepoint_compound_D4, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - 2D mCell D4 with compound'
) + theme(text = element_text(size = 15), 
          legend.text = element_text(size = 15), 
          aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_2D_mCell_timepoint_compound_D4 , file = 'GOpathways_2D_mCell_timepoint_compound_D4.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_2D_mCell_withcompound_D4toD0.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_2D_mCell_compound_D4 <- assay(rld_2D_mCell_timepoint_compound_D4)[rownames(assay(rld_2D_mCell_timepoint_compound_D4)) %in% 
                                                                          DEgenelist_2D_mCell_timepoint_compound_D4, c(1:4, 13:16)]
head(DE_normcounts_2D_mCell_compound_D4)
nrow(DE_normcounts_2D_mCell_compound_D4)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_2D_mCell_compound_D4 <- as.data.frame(t(scale(t(DE_normcounts_2D_mCell_compound_D4))))
head(scale_DE_normcounts_2D_mCell_compound_D4)


# color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_2D_mCell_compound_D4 <- Heatmap(scale_DE_normcounts_2D_mCell_compound_D4,
                                  col = color_ht,
                                  width = unit(10*ncol(scale_DE_normcounts_2D_mCell_compound_D4), 'mm'), 
                                  # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_compound_D4), 'mm'),
                                  name = 'Normalized counts', 
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE, 
                                  # show_row_names = FALSE,
                                  column_names_side = 'top', 
                                  # rect_gp = gpar(col = 'white', lwd = 1),
                                  # column_split = column_split_overlap,
                                  show_row_names = FALSE,
                                  row_names_max_width = max_text_width(
                                    rownames(scale_DE_normcounts_2D_mCell_compound_D4),
                                    gp = gpar(fontsize = 12))
)
ht_DE_2D_mCell_compound_D4




# Analysis of in vivo acute injury, find DE genes and pathways ------ 

# PCA plot of time points
library(ggplot2)
library(DESeq2)

rld_invivo_mCell_timepoint <- rlog(dds_invivo_mCell_timepoint)
pcaData_invivo <- plotPCA(rld_invivo_mCell_timepoint, intgroup = c('Timepoint'), returnData = T)
perVar_invivo <- round(100 * attr(pcaData_invivo, 'percentVar'))

head(pcaData_invivo)
# change the timepoint names, make it more professional
pcaData_invivo$group <- factor(c(rep('Healthy', 4), rep('Injury D1', 2), rep('Injury D3', 4), rep('Injury D7', 4)))
# adjust the order
pcaData_invivo$group <- factor(pcaData_invivo$group, levels = c('Healthy', 'Injury D1', 'Injury D3', 'Injury D7'))

ggplot(pcaData_invivo, aes(PC1, PC2, color = group)) + 
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", perVar_2D_solv[1], "% variance")) + 
  ylab(paste0("PC2: ", perVar_2D_solv[2], "% variance")) + 
  coord_fixed()+ 
  ggforce::geom_mark_ellipse(aes(fill = group)) + 
  theme(text = element_text(size = 20))  + ggtitle("PCA plot - in vivo acute injury")+
  xlim(-55, 50) +
  ylim(-25, 30) +
  
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 13/16
  )

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'PCA_in_vivo_acute_injury.tiff', 
  dpi = 300
)


## DE genes and GO pathways of in vivo D1 (acute) ===== 

summary(res_invivo_mCell_timepoint_acute)
  # out of 34595 with nonzero total read count
  # adjusted p-value < 0.05
  # LFC > 0 (up)       : 1175, 3.4%
  # LFC < 0 (down)     : 1587, 4.6%
  # outliers [1]       : 194, 0.56%
  # low counts [2]     : 11345, 33%
  # (mean count < 7)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results
res_invivo_mCell_timepoint_acute <- as.data.frame(res_invivo_mCell_timepoint_acute)
head(res_invivo_mCell_timepoint_acute)
sum(is.na(res_invivo_mCell_timepoint_acute$padj))
res_invivo_mCell_timepoint_acute <- na.omit(res_invivo_mCell_timepoint_acute)

DEgenelist_invivo_mCell_timepoint_acute <- rownames(res_invivo_mCell_timepoint_acute[res_invivo_mCell_timepoint_acute$padj < 0.05 & 
                                                                                   abs(res_invivo_mCell_timepoint_acute$log2FoldChange) >= 1, ])
head(DEgenelist_invivo_mCell_timepoint_acute)
length(DEgenelist_invivo_mCell_timepoint_acute)

DE_invivo_mCell_timepoint_acute <- res_invivo_mCell_timepoint_acute[res_invivo_mCell_timepoint_acute$padj < 0.05 & 
                                                                  abs(res_invivo_mCell_timepoint_acute$log2FoldChange) >=1, ]
head(DE_invivo_mCell_timepoint_acute)
nrow(DE_invivo_mCell_timepoint_acute)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
GOpathway_invivo_mCell_timepoint_acute <- enrichGO(gene = DEgenelist_invivo_mCell_timepoint_acute, 
                                                 OrgDb = org.Mm.eg.db, 
                                                 ont = 'BP', 
                                                 keyType = 'SYMBOL',
                                                 pAdjustMethod = 'BH', 
                                                 pvalueCutoff = 0.01,
                                                 qvalueCutoff = 0.05,
                                                 readable = TRUE)
head(GOpathway_invivo_mCell_timepoint_acute)
clusterProfiler::dotplot(GOpathway_invivo_mCell_timepoint_acute, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - mCell in vivo Treat acute injury D1') + theme(text = element_text(size = 15), 
                         legend.text = element_text(size = 15), 
                         aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_invivo_mCell_timepoint_acute, file = 'GOpathway_invivo_mCell_timepoint_acute.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_invivo_mCell_timepoint_acute.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_invivo_mCell_timepoint_acute <- assay(rld_invivo_mCell_timepoint_acute)[rownames(assay(rld_invivo_mCell_timepoint_acute)) %in% 
                                                                          DEgenelist_invivo_mCell_timepoint_acute, 1:6]
head(DE_normcounts_invivo_mCell_timepoint_acute)
nrow(DE_normcounts_invivo_mCell_timepoint_acute)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_invivo_mCell_timepoint_acute <- as.data.frame(t(scale(t(DE_normcounts_invivo_mCell_timepoint_acute))))
head(scale_DE_normcounts_invivo_mCell_timepoint_acute)


color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_invivo_mCell_timepoint_acute <- Heatmap(scale_DE_normcounts_invivo_mCell_timepoint_acute,
                                  col = color_ht,
                                  width = unit(10*ncol(scale_DE_normcounts_invivo_mCell_timepoint_acute), 'mm'), 
                                  # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_invivo_D1), 'mm'),
                                  name = 'Normalized counts', 
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE, 
                                  # show_row_names = FALSE,
                                  column_names_side = 'top', 
                                  # rect_gp = gpar(col = 'white', lwd = 1),
                                  # column_split = column_split_overlap,
                                  show_row_names = FALSE,
                                  row_names_max_width = max_text_width(
                                    rownames(scale_DE_normcounts_invivo_mCell_timepoint_acute),
                                    gp = gpar(fontsize = 12))
)
ht_DE_invivo_mCell_timepoint_acute


## DE genes and GO pathways of in vivo D3 ===== 

summary(res_invivo_mCell_timepoint_D3)
# out of 34595 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2319, 6.7%
# LFC < 0 (down)     : 1353, 3.9%
# outliers [1]       : 194, 0.56%
# low counts [2]     : 10680, 31%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
res_invivo_mCell_timepoint_D3 <- as.data.frame(res_invivo_mCell_timepoint_D3)
head(res_invivo_mCell_timepoint_D3)
sum(is.na(res_invivo_mCell_timepoint_D3$padj))
res_invivo_mCell_timepoint_D3 <- na.omit(res_invivo_mCell_timepoint_D3)

DEgenelist_invivo_mCell_timepoint_D3 <- rownames(res_invivo_mCell_timepoint_D3[res_invivo_mCell_timepoint_D3$padj < 0.05 & 
                                                                                       abs(res_invivo_mCell_timepoint_D3$log2FoldChange) >= 1, ])
head(DEgenelist_invivo_mCell_timepoint_D3)
length(DEgenelist_invivo_mCell_timepoint_D3)

DE_invivo_mCell_timepoint_D3 <- res_invivo_mCell_timepoint_D3[res_invivo_mCell_timepoint_D3$padj < 0.05 & 
                                                                      abs(res_invivo_mCell_timepoint_D3$log2FoldChange) >=1, ]
head(DE_invivo_mCell_timepoint_D3)
nrow(DE_invivo_mCell_timepoint_D3)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
GOpathway_invivo_mCell_timepoint_D3 <- enrichGO(gene = DEgenelist_invivo_mCell_timepoint_D3, 
                                                   OrgDb = org.Mm.eg.db, 
                                                   ont = 'BP', 
                                                   keyType = 'SYMBOL',
                                                   pAdjustMethod = 'BH', 
                                                   pvalueCutoff = 0.01,
                                                   qvalueCutoff = 0.05,
                                                   readable = TRUE)
head(GOpathway_invivo_mCell_timepoint_D3)
clusterProfiler::dotplot(GOpathway_invivo_mCell_timepoint_D3, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - mCell in vivo Treat acute injury D3') + theme(text = element_text(size = 15), 
                                                                                             legend.text = element_text(size = 15), 
                                                                                             aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_invivo_mCell_timepoint_D3, file = 'GOpathway_invivo_mCell_timepoint_D3.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_invivo_mCell_timepoint_D3.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_invivo_mCell_timepoint_D3 <- assay(rld_invivo_mCell_timepoint_D3)[rownames(assay(rld_invivo_mCell_timepoint_D3)) %in% 
                                                                                        DEgenelist_invivo_mCell_timepoint_D3, c(1:4, 7:10)]
head(DE_normcounts_invivo_mCell_timepoint_D3)
nrow(DE_normcounts_invivo_mCell_timepoint_D3)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_invivo_mCell_timepoint_D3 <- as.data.frame(t(scale(t(DE_normcounts_invivo_mCell_timepoint_D3))))
head(scale_DE_normcounts_invivo_mCell_timepoint_D3)


color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_invivo_mCell_timepoint_D3 <- Heatmap(scale_DE_normcounts_invivo_mCell_timepoint_D3,
                                              col = color_ht,
                                              width = unit(10*ncol(scale_DE_normcounts_invivo_mCell_timepoint_D3), 'mm'), 
                                              # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_invivo_D1), 'mm'),
                                              name = 'Normalized counts', 
                                              cluster_columns = FALSE,
                                              cluster_rows = TRUE, 
                                              # show_row_names = FALSE,
                                              column_names_side = 'top', 
                                              # rect_gp = gpar(col = 'white', lwd = 1),
                                              # column_split = column_split_overlap,
                                              show_row_names = FALSE,
                                              row_names_max_width = max_text_width(
                                                rownames(scale_DE_normcounts_invivo_mCell_timepoint_D3),
                                                gp = gpar(fontsize = 12))
)
ht_DE_invivo_mCell_timepoint_D3

## DE genes and GO pathways of in vivo D7 (recover) =====

summary(res_invivo_mCell_timepoint_recover)
  # out of 34595 with nonzero total read count
  # adjusted p-value < 0.05
  # LFC > 0 (up)       : 1683, 4.9%
  # LFC < 0 (down)     : 960, 2.8%
  # outliers [1]       : 194, 0.56%
  # low counts [2]     : 10680, 31%
  # (mean count < 6)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results
res_invivo_mCell_timepoint_recover <- as.data.frame(res_invivo_mCell_timepoint_recover)
head(res_invivo_mCell_timepoint_recover)
sum(is.na(res_invivo_mCell_timepoint_recover$padj))
res_invivo_mCell_timepoint_recover <- na.omit(res_invivo_mCell_timepoint_recover)

DEgenelist_invivo_mCell_timepoint_recover <- rownames(res_invivo_mCell_timepoint_recover[res_invivo_mCell_timepoint_recover$padj < 0.05 & 
                                                                                       abs(res_invivo_mCell_timepoint_recover$log2FoldChange) >= 1, ])
head(DEgenelist_invivo_mCell_timepoint_recover)
length(DEgenelist_invivo_mCell_timepoint_recover)

DE_invivo_mCell_timepoint_recover <- res_invivo_mCell_timepoint_recover[res_invivo_mCell_timepoint_recover$padj < 0.05 & 
                                                                      abs(res_invivo_mCell_timepoint_recover$log2FoldChange) >=1, ]
head(DE_invivo_mCell_timepoint_recover)
nrow(DE_invivo_mCell_timepoint_recover)

# GO pathways
library(clusterProfiler)
library(org.Mm.eg.db)
GOpathway_invivo_mCell_timepoint_recover <- enrichGO(gene = DEgenelist_invivo_mCell_timepoint_recover, 
                                                   OrgDb = org.Mm.eg.db, 
                                                   ont = 'BP', 
                                                   keyType = 'SYMBOL',
                                                   pAdjustMethod = 'BH', 
                                                   pvalueCutoff = 1,
                                                   qvalueCutoff = 1,
                                                   readable = TRUE)
head(GOpathway_invivo_mCell_timepoint_recover)
clusterProfiler::dotplot(GOpathway_invivo_mCell_timepoint_recover, showCategory = 20, 
                         font.size = 14, 
                         label_format = 80, 
                         title = 'GO analysis - mCell in vivo Treat acute injury D7 - recover') + theme(text = element_text(size = 15), 
                                                                                             legend.text = element_text(size = 15), 
                                                                                             aspect.ratio = 2/1) 

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
library(xlsx)
write.xlsx(GOpathway_invivo_mCell_timepoint_recover, file = 'GOpathway_invivo_mCell_timepoint_recover.xlsx', 
           row.names = TRUE, col.names = TRUE)


# ERK1 and ERK2: mediate Cell phenotypic alternation during chronic liver injury (10.1016/j.ymthe.2018.08.016)

# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'GOpathway_invivo_mCell_timepoint_recover.tiff', 
  plot = last_plot(),
  dpi = 300, 
  width = 14, 
  height = 7
)

# Heatmap of DE genes
library(ComplexHeatmap)
library(circlize)
DE_normcounts_invivo_mCell_timepoint_recover <- assay(rld_invivo_mCell_timepoint_recover)[rownames(assay(rld_invivo_mCell_timepoint_recover)) %in% 
                                                                                        DEgenelist_invivo_mCell_timepoint_recover, c(1:4, 11:14)]
head(DE_normcounts_invivo_mCell_timepoint_recover)
nrow(DE_normcounts_invivo_mCell_timepoint_recover)

# scale and calculate mean of the DE norm counts
scale_DE_normcounts_invivo_mCell_timepoint_recover <- as.data.frame(t(scale(t(DE_normcounts_invivo_mCell_timepoint_recover))))
head(scale_DE_normcounts_invivo_mCell_timepoint_recover)


color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

ht_DE_invivo_mCell_timepoint_recover <- Heatmap(scale_DE_normcounts_invivo_mCell_timepoint_recover,
                                              col = color_ht,
                                              width = unit(10*ncol(scale_DE_normcounts_invivo_mCell_timepoint_recover), 'mm'), 
                                              # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_invivo_D1), 'mm'),
                                              name = 'Normalized counts', 
                                              cluster_columns = FALSE,
                                              cluster_rows = TRUE, 
                                              # show_row_names = FALSE,
                                              column_names_side = 'top', 
                                              # rect_gp = gpar(col = 'white', lwd = 1),
                                              # column_split = column_split_overlap,
                                              show_row_names = FALSE,
                                              row_names_max_width = max_text_width(
                                                rownames(scale_DE_normcounts_invivo_mCell_timepoint_recover),
                                                gp = gpar(fontsize = 12))
)
ht_DE_invivo_mCell_timepoint_recover

# venn diagram of DE genes in vivo acute injury
# library(ggvenn)
# DEgeneset_invivo <- list(`Acute injury` = DEgenelist_invivo_mCell_timepoint_acute, 
#                          `Recovery` = DEgenelist_invivo_mCell_timepoint_recover)
# 
# ggvenn(DEgeneset_invivo, c('Acute injury', 'Recovery'))

# Compare the change of acute injury and after seeding -------
library(ggvenn)
DEgeneset_invivo_2D <- list(`Acute injury` = DEgenelist_invivo_mCell_timepoint_acute,
                            `2D culture D1 without compound` = DEgenelist_2D_mCell_timepoint_solv_D1, 
                            `2D culture D1 with compound` = DEgenelist_2D_mCell_timepoint_compound_D1, 
                            `2D culture D4 without compound` = DEgenelist_2D_mCell_timepoint_solv_D4, 
                            `2D culture D4 with compound` = DEgenelist_2D_mCell_timepoint_compound_D4, 
                            `Solv D2` = DEgenelist_2D_mCell_timepoint_solv_D2, 
                            `compound D2` = DEgenelist_2D_mCell_timepoint_compound_D2
                            )
ggvenn(DEgeneset_invivo_2D, c('Acute injury', '2D culture D1 without compound'))
# save high resolution images
setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'venn_acuteinjury_D1solv.tiff', 
  dpi = 300
)
ggvenn(DEgeneset_invivo_2D, c('Acute injury', '2D culture D1 with compound'))
ggvenn(DEgeneset_invivo_2D, c('2D culture D1 with compound', '2D culture D1 without compound'))
ggvenn(DEgeneset_invivo_2D, c('Acute injury', '2D culture D4 without compound'))  
ggvenn(DEgeneset_invivo_2D, c('Solv D2', 'compound D2'))
       
# <!Interesting part!> Compare the pathways in different conditions ------------

# list_upset_pathway <- list(`2D D1 withoutcompound` = GOpathway_2D_mCell_timepoint_solv_D1[1:20, ]$Description, 
#                            `2D D2 withoutcompound` = GOpathway_2D_mCell_timepoint_solv_D2[1:20, ]$Description, 
#                            `2D D4 withoutcompound` = GOpathway_2D_mCell_timepoint_solv_D4[1:20, ]$Description, 
#                            `2D D1 withcompound` = GOpathway_2D_mCell_timepoint_compound_D1[1:20, ]$Description,
#                            `2D D2 withcompound` = GOpathway_2D_mCell_timepoint_compound_D2[1:20, ]$Description,
#                            `2D D4 withcompound` = GOpathway_2D_mCell_timepoint_compound_D4[1:20, ]$Description,
#                            `Acute injury D1` = GOpathway_invivo_mCell_timepoint_acute[1:20, ]$Description, 
#                            `Acute injury D3` = GOpathway_invivo_mCell_timepoint_D3[1:20, ]$Description, 
#                            `Acute injury D7` = GOpathway_invivo_mCell_timepoint_recover[1:20, ]$Description
#                            )

# ggvenn and upset plot cannot completely see the overlap (in upset plot they are separated into different groups)
# Directly subset the overlap and draw plot (maybe heatmap?)
top_GOpathway_invivo_mCell_timepoint_acute <- GOpathway_invivo_mCell_timepoint_acute[1:20, ]
top_GOpathway_invivo_mCell_timepoint_D3 <- GOpathway_invivo_mCell_timepoint_D3[1:20, ]
top_GOpathway_invivo_mCell_timepoint_recover <- GOpathway_invivo_mCell_timepoint_recover[1:20, ]

top_GOpathway_2D_mCell_timepoint_solv_D1 <- GOpathway_2D_mCell_timepoint_solv_D1[1:20, ]
top_GOpathway_2D_mCell_timepoint_solv_D2 <- GOpathway_2D_mCell_timepoint_solv_D2[1:20, ]
top_GOpathway_2D_mCell_timepoint_solv_D4 <- GOpathway_2D_mCell_timepoint_solv_D4[1:20, ]

top_GOpathway_2D_mCell_timepoint_compound_D1 <- GOpathway_2D_mCell_timepoint_compound_D1[1:20, ]
top_GOpathway_2D_mCell_timepoint_compound_D2 <- GOpathway_2D_mCell_timepoint_compound_D2[1:20, ]
top_GOpathway_2D_mCell_timepoint_compound_D4 <- GOpathway_2D_mCell_timepoint_compound_D4[1:20, ]

# rm(GOpathway_invivo_mCell_timepoint_acute)
# rm(GOpathway_invivo_mCell_timepoint_D3)
# rm(GOpathway_invivo_mCell_timepoint_recover)
# 
# rm(GOpathway_2D_mCell_timepoint_solv_D1)
# rm(GOpathway_2D_mCell_timepoint_solv_D2)
# rm(GOpathway_2D_mCell_timepoint_solv_D4)
# 
# rm(GOpathway_2D_mCell_timepoint_compound_D1)
# rm(GOpathway_2D_mCell_timepoint_compound_D2)
# rm(GOpathway_2D_mCell_timepoint_compound_D4)

## acute in vivo compared to 2D groups =========
library(dplyr)
overlap_pathways_2D_D1_Acute <- inner_join(top_GOpathway_2D_mCell_timepoint_solv_D1, top_GOpathway_invivo_mCell_timepoint_acute, 
                                           by = 'ID')
head(overlap_pathways_2D_D1_Acute)

overlap_pathways_2D_D2_Acute <- top_GOpathway_invivo_mCell_timepoint_acute[
  top_GOpathway_invivo_mCell_timepoint_acute$Description %in% top_GOpathway_2D_mCell_timepoint_solv_D2$Description, 
]
head(overlap_pathways_2D_D2_Acute)

overlap_pathways_2D_D4_Acute <- top_GOpathway_invivo_mCell_timepoint_acute[
  top_GOpathway_invivo_mCell_timepoint_acute$Description %in% top_GOpathway_2D_mCell_timepoint_solv_D4$Description, 
]
head(overlap_pathways_2D_D4_Acute)
# there is no overlapped pathway between 2D culture D4 and in vivo acute injury (D1)

overlap_pathways_2D_compound_D1_Acute <- top_GOpathway_invivo_mCell_timepoint_acute[
  top_GOpathway_invivo_mCell_timepoint_acute$Description %in% top_GOpathway_2D_mCell_timepoint_compound_D1$Description, 
]
head(overlap_pathways_2D_compound_D1_Acute)

overlap_pathways_2D_compound_D2_Acute <- top_GOpathway_invivo_mCell_timepoint_acute[
  top_GOpathway_invivo_mCell_timepoint_acute$Description %in% top_GOpathway_2D_mCell_timepoint_compound_D2$Description, 
]
head(overlap_pathways_2D_compound_D2_Acute)

overlap_pathways_2D_compound_D4_Acute <- top_GOpathway_invivo_mCell_timepoint_acute[
  top_GOpathway_invivo_mCell_timepoint_acute$Description %in% top_GOpathway_2D_mCell_timepoint_compound_D4$Description, 
]
head(overlap_pathways_2D_compound_D4_Acute)
# there is no overlapped pathway between 2D culture D4 and in vivo acute injury (D1)

# Check the genes






print('123')






# Draw the overlap gene heatmap of acute injury (D1) and 2D culture D1 without compound ======
library(ComplexHeatmap)
library(circlize)
Overlapgene_acute_2Dsolv_D1 <- DE_invivo_mCell_timepoint_acute[rownames(DE_invivo_mCell_timepoint_acute) %in% 
                                  rownames(DE_2D_mCell_timepoint_solv_D1),]
Overlapgene_acute_2Dsolv_2D <- assay(rld_2D_mCell_solv)[rownames(assay(rld_2D_mCell_solv)) %in% rownames(Overlapgene_acute_2Dsolv_D1), 
                                                        1:8]
Overlapgene_acute_2Dsolv_invivo <- assay(rld_invivo_mCell_timepoint)[rownames(assay(rld_invivo_mCell_timepoint)) %in% 
                                                                rownames(Overlapgene_acute_2Dsolv), 1:6]
head(Overlapgene_acute_2Dsolv_2D)
head(Overlapgene_acute_2Dsolv_invivo)

color_ht <- colorRamp2(c(-1, 0, 1), c("#1f78b4", "white", "#ff7f00"))

# scale and calculate mean of the DE norm counts
scale_Overlapgene_acute_2Dsolv_2D <- as.data.frame(t(scale(t(Overlapgene_acute_2Dsolv_2D))))
scale_Overlapgene_acute_2Dsolv_invivo <- as.data.frame(t(scale(t(Overlapgene_acute_2Dsolv_invivo))))
head(scale_Overlapgene_acute_2Dsolv_2D)
head(scale_Overlapgene_acute_2Dsolv_invivo)

library(dplyr)

ht_Overlapgene_acute_2Dsolv_2D <- Heatmap(scale_Overlapgene_acute_2Dsolv_2D,
                                          col = color_ht,
                                          width = unit(10*ncol(scale_Overlapgene_acute_2Dsolv_2D), 'mm'), 
                                          # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_solv_D1), 'mm'),
                                          # name = 'Normalized counts', 
                                          cluster_columns = FALSE,
                                          cluster_rows = TRUE, 
                                          # show_row_names = FALSE,
                                          column_names_side = 'top', 
                                          # rect_gp = gpar(col = 'white', lwd = 1),
                                          # column_split = column_split_overlap,
                                          show_row_names = FALSE,
                                          row_names_max_width = max_text_width(
                                            rownames(scale_Overlapgene_acute_2Dsolv_2D),
                                            gp = gpar(fontsize = 12)))
ht_Overlapgene_acute_2Dsolv_2D


ht_Overlapgene_acute_2Dsolv_invivo <- Heatmap(scale_Overlapgene_acute_2Dsolv_invivo,
                                          col = color_ht,
                                          width = unit(10*ncol(scale_Overlapgene_acute_2Dsolv_invivo), 'mm'), 
                                          # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_solv_D1), 'mm'),
                                          # name = 'Normalized counts', 
                                          cluster_columns = FALSE,
                                          cluster_rows = TRUE, 
                                          # show_row_names = FALSE,
                                          column_names_side = 'top', 
                                          # rect_gp = gpar(col = 'white', lwd = 1),
                                          # column_split = column_split_overlap,
                                          show_row_names = FALSE,
                                          row_names_max_width = max_text_width(
                                            rownames(scale_Overlapgene_acute_2Dsolv_invivo),
                                            gp = gpar(fontsize = 12)))
ht_Overlapgene_acute_2Dsolv_invivo

# Check the overlap genes with the same trend
Overlapgene_acute_2Dsolv_up <- Overlapgene_acute_2Dsolv[Overlapgene_acute_2Dsolv$log2FoldChange>0,]
Overlapgene_acute_2Dsolv_down <- Overlapgene_acute_2Dsolv[Overlapgene_acute_2Dsolv$log2FoldChange<0,]
head(Overlapgene_acute_2Dsolv_up)
nrow(Overlapgene_acute_2Dsolv_up)
nrow(Overlapgene_acute_2Dsolv_down)

# up
Overlapgene_acute_2Dsolv_up_counts_2D <- Overlapgene_acute_2Dsolv_2D[rownames(Overlapgene_acute_2Dsolv_2D) %in% 
                                                                    rownames(Overlapgene_acute_2Dsolv_up), ]
head(Overlapgene_acute_2Dsolv_up_counts_2D)
nrow(Overlapgene_acute_2Dsolv_up_counts_2D)
Overlapgene_acute_2Dsolv_up_counts_invivo <- Overlapgene_acute_2Dsolv_invivo[rownames(Overlapgene_acute_2Dsolv_invivo) %in% 
                                                                               rownames(Overlapgene_acute_2Dsolv_up), ]
head(Overlapgene_acute_2Dsolv_up_counts_invivo)
nrow(Overlapgene_acute_2Dsolv_up_counts_invivo)

Overlapgene_acute_2Dsolv_up_counts <- cbind(Overlapgene_acute_2Dsolv_up_counts_2D, Overlapgene_acute_2Dsolv_up_counts_invivo)
head(Overlapgene_acute_2Dsolv_up_counts)
nrow(Overlapgene_acute_2Dsolv_up_counts)
# rm(Overlapgene_acute_2Dsolv_up_counts_2D)
# rm(Overlapgene_acute_2Dsolv_up_counts_invivo)

#down
Overlapgene_acute_2Dsolv_down_counts_2D <- Overlapgene_acute_2Dsolv_2D[rownames(Overlapgene_acute_2Dsolv_2D) %in% 
                                                                       rownames(Overlapgene_acute_2Dsolv_down), ]
head(Overlapgene_acute_2Dsolv_down_counts_2D)
nrow(Overlapgene_acute_2Dsolv_down_counts_2D)
Overlapgene_acute_2Dsolv_down_counts_invivo <- Overlapgene_acute_2Dsolv_invivo[rownames(Overlapgene_acute_2Dsolv_invivo) %in% 
                                                                               rownames(Overlapgene_acute_2Dsolv_down), ]
head(Overlapgene_acute_2Dsolv_down_counts_invivo)
nrow(Overlapgene_acute_2Dsolv_down_counts_invivo)

Overlapgene_acute_2Dsolv_down_counts <- cbind(Overlapgene_acute_2Dsolv_down_counts_2D, Overlapgene_acute_2Dsolv_down_counts_invivo)
head(Overlapgene_acute_2Dsolv_down_counts)
nrow(Overlapgene_acute_2Dsolv_down_counts)

# rm(Overlapgene_acute_2Dsolv_down_counts_2D)
# rm(Overlapgene_acute_2Dsolv_down_counts_invivo)


# Heatmap of the same trend genes
scale_Overlapgene_acute_2Dsolv_up_counts <- as.data.frame(t(scale(t(Overlapgene_acute_2Dsolv_up_counts))))
head(scale_Overlapgene_acute_2Dsolv_up_counts)
scale_Overlapgene_acute_2Dsolv_up_counts_mean <- mutate(scale_Overlapgene_acute_2Dsolv_up_counts, 
                                                        `Fresh isolated` = apply(scale_Overlapgene_acute_2Dsolv_up_counts[, 1:4], 1, mean))
scale_Overlapgene_acute_2Dsolv_up_counts_mean <- mutate(scale_Overlapgene_acute_2Dsolv_up_counts_mean, 
                                                        `Day 1` = apply(scale_Overlapgene_acute_2Dsolv_up_counts[, 5:8], 1, mean))
scale_Overlapgene_acute_2Dsolv_up_counts_mean <- mutate(scale_Overlapgene_acute_2Dsolv_up_counts_mean, 
                                                        `Healthy` = apply(scale_Overlapgene_acute_2Dsolv_up_counts[, 9:12], 1, mean))
scale_Overlapgene_acute_2Dsolv_up_counts_mean <- mutate(scale_Overlapgene_acute_2Dsolv_up_counts_mean, 
                                                        `Acute injury - D1` = apply(scale_Overlapgene_acute_2Dsolv_up_counts[, 13:14], 1, mean))
scale_Overlapgene_acute_2Dsolv_up_counts_mean <- scale_Overlapgene_acute_2Dsolv_up_counts_mean[ , c('Fresh isolated', 
                                                                                                    'Day 1', 'Healthy', 'Acute injury - D1')]
head(scale_Overlapgene_acute_2Dsolv_up_counts_mean)

Overlapgene_acute_2Dsolv_up  <- Heatmap(scale_Overlapgene_acute_2Dsolv_up_counts_mean,
                                              col = color_ht,
                                              width = unit(10*ncol(scale_Overlapgene_acute_2Dsolv_up_counts_mean), 'mm'), 
                                              # height = unit(5* nrow(scale_DE_normcounts_2D_mCell_solv_D1), 'mm'),
                                              # name = 'Normalized counts', 
                                              cluster_columns = FALSE,
                                              cluster_rows = TRUE, 
                                              # show_row_names = FALSE,
                                              column_names_side = 'top', 
                                              # rect_gp = gpar(col = 'white', lwd = 1),
                                              # column_split = column_split_overlap,
                                              show_row_names = FALSE,
                                              row_names_max_width = max_text_width(
                                                rownames(scale_Overlapgene_acute_2Dsolv_up_counts_mean),
                                                gp = gpar(fontsize = 12)))
Overlapgene_acute_2Dsolv_up





# Specific gene markers -------------
rld_counts_invitro <- assay(rld_2D_mCell_timepoint)
rld_counts_invivo <- assay(rld_invivo_mCell_timepoint)
rld_counts_3D <- assay(rld_3D)

Typical_Cell_genes <- c('Stab2', 'Clec4g', 'Lyve1', 'Clec1b', 'Cd36', 'Pecam1','Fcgr2b', 'Flt4')

Typical_Cell_gene_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                      Typical_Cell_genes, ]
Typical_Cell_gene_rld_counts_2D_mCell_invitro <- as.data.frame(t(Typical_Cell_gene_rld_counts_2D_mCell_invitro))
Typical_Cell_gene_rld_counts_2D_mCell_invitro$group <- c(rep('D0', 4), rep('D1_Solv', 4), 
                                                         rep('D2_Solv', 4), rep('D4_Solv',4), 
                                                         rep('D1_compound', 4), 
                                                         rep('D2_compound', 4), rep('D4_compound',4)
                                                         )
head(Typical_Cell_gene_rld_counts_2D_mCell_invitro)

Typical_Cell_gene_rld_counts_2D_mCell_solv <- Typical_Cell_gene_rld_counts_2D_mCell_invitro[1:16,]
head(Typical_Cell_gene_rld_counts_2D_mCell_solv)

Typical_Cell_gene_rld_counts_2D_mCell_compound <- Typical_Cell_gene_rld_counts_2D_mCell_invitro[c(1:4, 17:28), ]
head(Typical_Cell_gene_rld_counts_2D_mCell_compound)

Typical_Cell_gene_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                           Typical_Cell_genes, ]
Typical_Cell_gene_rld_counts_invivo <- as.data.frame(t(Typical_Cell_gene_rld_counts_invivo))
Typical_Cell_gene_rld_counts_invivo$group <- c(rep('0', 4), rep('1', 2), rep('3', 4), rep('7',4))
head(Typical_Cell_gene_rld_counts_invivo)


Typical_Cell_gene_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                           Typical_Cell_genes, ]
Typical_Cell_gene_rld_counts_3D <- as.data.frame(t(Typical_Cell_gene_rld_counts_3D))
Typical_Cell_gene_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Typical_Cell_gene_rld_counts_3D)


library(ggplot2)
library(reshape2)

# invitro Cell Typical gene markers
Typical_Cell_gene_rld_counts_2D_mCell_solv <- melt(Typical_Cell_gene_rld_counts_2D_mCell_solv, id.vars = 'group')
Typical_Cell_gene_rld_counts_2D_mCell_compound <- melt(Typical_Cell_gene_rld_counts_2D_mCell_compound, id.vars = 'group')
head(Typical_Cell_gene_rld_counts_2D_mCell_solv)
tail(Typical_Cell_gene_rld_counts_2D_mCell_solv)
tail(Typical_Cell_gene_rld_counts_2D_mCell_compound)
Typical_Cell_gene_rld_counts_2D_mCell_solv$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
                                                           levels = c('0', '1', '2', '4'))
Typical_Cell_gene_rld_counts_2D_mCell_compound$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
                                                           levels = c('0', '1', '2', '4'))
colnames(Typical_Cell_gene_rld_counts_2D_mCell_solv) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
colnames(Typical_Cell_gene_rld_counts_2D_mCell_compound) <- c('Timepoint', 'variable', 'log2(normalized_counts)')

ggplot(data = Typical_Cell_gene_rld_counts_2D_mCell_solv, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
  geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Typical Cell gene markers - 2D without compound') + 
  geom_jitter(shape = 16, position = position_jitter(0.2))

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'boxplot_2D_solv_TypicalCellmarkers.tiff', 
  dpi = 300
)

ggplot(data = Typical_Cell_gene_rld_counts_2D_mCell_compound, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
  geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Typical Cell gene markers - 2D with compound')+ 
  geom_jitter(shape = 16, position = position_jitter(0.2))

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'boxplot_2D_compound_TypicalCellmarkers.tiff', 
  dpi = 300
)



# invivo Cell Typical gene markers
Typical_Cell_gene_rld_counts_invivo <- melt(Typical_Cell_gene_rld_counts_invivo, id.vars = 'group')
Typical_Cell_gene_rld_counts_invivo$group <- factor(Typical_Cell_gene_rld_counts_invivo$group, 
                                                      levels = c('0', '1', '3', '7'))
colnames(Typical_Cell_gene_rld_counts_invivo) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
ggplot(data = Typical_Cell_gene_rld_counts_invivo, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
  geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Typical Cell gene markers - in vivo acute injury')+ 
  geom_jitter(shape = 16, position = position_jitter(0.2))

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'boxplot_invivo_acute_injury_TypicalCellmarkers.tiff', 
  dpi = 300
)


# 3D Cell Typical gene markers
Typical_Cell_gene_rld_counts_3D <- melt(Typical_Cell_gene_rld_counts_3D, id.vars = 'group')
Typical_Cell_gene_rld_counts_3D$group <- factor(Typical_Cell_gene_rld_counts_3D$group, 
                                                               levels = c('0', '4', '7', '10', '14'))
colnames(Typical_Cell_gene_rld_counts_3D) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
ggplot(data = Typical_Cell_gene_rld_counts_3D, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
  geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Typical Cell gene markers - 3D')+ 
  geom_jitter(shape = 16, position = position_jitter(0.2))

setwd('D:/foldername/Cell_Dedif/20230906/plots/')
ggplot2::ggsave(
  filename = 'boxplot_3D_TypicalCellmarkers.tiff', 
  dpi = 300
)






# Cell signatures from paper-------------------
Gene_Cell_signature_D7high <- c('Hspa1a', 'Hspa1b', 'Ecm1', 'Sema6a', 'Gm49369', 
                                'Timp2', 'Hipk2', 'Reln', 'Il33', 'Ctsd','Fcgr3', 
                                'Nr2f1', 'Pltp', 'Acp5')
Gene_Cell_signature_D1high <- c('Ccl6', 'Cd14', 'Ccl9', 'Insr', 'Adm', 'Nid1', 
                                'Timp1', 'Akap12', 'Lyve1', 'Clec1b', 'F2r', 
                                'Tfpi2', 'Cd4', 'Il6st', 'Crhbp', 'Cd36', 
                                'Stab1')
Gene_Cell_signature_Healthyhigh <-c('Clec4g','Stab2','Oit3','Gpr182','Fcgr2b','Mrc1',
                                    'Kdr','Phactr2','Nrp1','Maf','Lifr','Snx5','Hibch',
                                    'Cd209g','Cd209f','Dab2','Ntn4','Plpp3','Npl') 
Gene_Cell_damage_signature <- c('Fabp5','Vwf','Vwa1','Fabp4')




# Gene signatures D7 high - recovery signature =====
Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                      Gene_Cell_signature_D7high, ]
# Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro <- as.data.frame(t(Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro))
# Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro$group <- c(rep('D0', 4), rep('D1_Solv', 4), 
#                                                                    rep('D2_Solv', 4), rep('D4_Solv',4), 
#                                                                    rep('D1_compound', 4), 
#                                                                    rep('D2_compound', 4), rep('D4_compound',4)
#                                                                   )
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro)

Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv <- Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv)

Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound <- Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound)


Gene_Cell_signature_D7high_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                    Gene_Cell_signature_D7high, ]
# Gene_Cell_signature_D7high_rld_counts_invivo <- as.data.frame(t(Gene_Cell_signature_D7high_rld_counts_invivo))
# Gene_Cell_signature_D7high_rld_counts_invivo$group <- c(rep('0', 4), rep('1', 2), rep('3', 4), rep('7',4))
head(Gene_Cell_signature_D7high_rld_counts_invivo)


Gene_Cell_signature_D7high_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                        Gene_Cell_signature_D7high, ]
# Gene_Cell_signature_D7high_rld_counts_3D <- as.data.frame(t(Gene_Cell_signature_D7high_rld_counts_3D))
# Gene_Cell_signature_D7high_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_Cell_signature_D7high_rld_counts_3D)


library(ggplot2)
library(reshape2)

# invitro Gene signatures D7 high

Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv))))
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv)

Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv, 
                         `D0` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean, 
                         `D1` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean, 
                         `D2` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean, 
                         `D4` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean <- Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean)



Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound))))
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound)

Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound, 
                                                                   `D0` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean, 
                                                                   `D1` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean, 
                                                                   `D2` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean, 
                                                                   `D4` = apply(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean <- Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean)



Gene_Cell_signature_D7high_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_Cell_signature_D7high_rld_counts_invivo))))
head(Gene_Cell_signature_D7high_rld_counts_invivo)

Gene_Cell_signature_D7high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_invivo, 
                                                                   `D0` = apply(Gene_Cell_signature_D7high_rld_counts_invivo[1:4], 1, mean))
Gene_Cell_signature_D7high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_invivo_mean, 
                                                                   `D1` = apply(Gene_Cell_signature_D7high_rld_counts_invivo[5: 6], 1, mean))
Gene_Cell_signature_D7high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_invivo_mean, 
                                                                   `D3` = apply(Gene_Cell_signature_D7high_rld_counts_invivo[7: 10], 1, mean))
Gene_Cell_signature_D7high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_invivo_mean, 
                                                                   `D7` = apply(Gene_Cell_signature_D7high_rld_counts_invivo[11: 14], 1, mean))

Gene_Cell_signature_D7high_rld_counts_invivo_mean <- Gene_Cell_signature_D7high_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_Cell_signature_D7high_rld_counts_invivo_mean)


Gene_Cell_signature_D7high_rld_counts_3D <- as.data.frame(t(scale(t(Gene_Cell_signature_D7high_rld_counts_3D))))

Gene_Cell_signature_D7high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_3D, 
                        `D0` = apply(Gene_Cell_signature_D7high_rld_counts_3D[1: 3], 1, mean))
Gene_Cell_signature_D7high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_3D_mean, 
                        `D4` = apply(Gene_Cell_signature_D7high_rld_counts_3D_mean[4: 6], 1, mean))
Gene_Cell_signature_D7high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_3D_mean, 
                        `D7` = apply(Gene_Cell_signature_D7high_rld_counts_3D_mean[7: 9], 1, mean))
Gene_Cell_signature_D7high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_3D_mean, 
                        `D10` = apply(Gene_Cell_signature_D7high_rld_counts_3D_mean[10: 12], 1, mean))
Gene_Cell_signature_D7high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D7high_rld_counts_3D_mean, 
                        `D14` = apply(Gene_Cell_signature_D7high_rld_counts_3D_mean[13: 15], 1, mean))
Gene_Cell_signature_D7high_rld_counts_3D_mean <- Gene_Cell_signature_D7high_rld_counts_3D_mean[,16: 20]
head(Gene_Cell_signature_D7high_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_D7high_solv <- Heatmap(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean,
                       col = color_ht,
                       width = unit(8*ncol(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean), 'mm'), 
                       height = unit(6* nrow(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean), 'mm'),
                       row_title = 'Recovery gene signature', 
                       cluster_columns = FALSE,
                       cluster_rows = FALSE, 
                       # cluster_rows = FALSE,
                       # show_row_names = FALSE,
                       column_names_side = 'top', 
                       rect_gp = gpar(col = 'white', lwd = 1), 
                       show_heatmap_legend = FALSE
                       # column_split = column_split_overlap,
                       # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv_mean,
                       #                                      gp = gpar(fontsize = 12))    )
)

ht_D7high_compound <- Heatmap(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean,
                          col = color_ht,
                          width = unit(8*ncol(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean), 'mm'), 
                          height = unit(6* nrow(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean), 'mm'),
                          row_title = 'Recovery gene signature', 
                          cluster_columns = FALSE,
                          cluster_rows = FALSE, 
                          # cluster_rows = FALSE,
                          # show_row_names = FALSE,
                          column_names_side = 'top', 
                          rect_gp = gpar(col = 'white', lwd = 1), 
                          show_heatmap_legend = FALSE
                          # column_split = column_split_overlap,
                          # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound_mean,
                          #                                      gp = gpar(fontsize = 12))    )
)

ht_D7high_invivo <- Heatmap(Gene_Cell_signature_D7high_rld_counts_invivo_mean,
                          col = color_ht,
                          width = unit(8*ncol(Gene_Cell_signature_D7high_rld_counts_invivo_mean), 'mm'), 
                          height = unit(6* nrow(Gene_Cell_signature_D7high_rld_counts_invivo_mean), 'mm'),
                          row_title = 'Recovery gene signature', 
                          cluster_columns = FALSE,
                          cluster_rows = FALSE, 
                          # cluster_rows = FALSE,
                          # show_row_names = FALSE,
                          column_names_side = 'top', 
                          rect_gp = gpar(col = 'white', lwd = 1), 
                          show_heatmap_legend = FALSE
                          # column_split = column_split_overlap,
                          # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D7high_rld_counts_2D_mCell_invivo_mean,
                          #                                      gp = gpar(fontsize = 12))    )
)

ht_D7high_3D <- Heatmap(Gene_Cell_signature_D7high_rld_counts_3D_mean,
                            col = color_ht,
                            width = unit(8*ncol(Gene_Cell_signature_D7high_rld_counts_3D_mean), 'mm'), 
                            height = unit(6* nrow(Gene_Cell_signature_D7high_rld_counts_3D_mean), 'mm'),
                            row_title = 'Recovery gene signature', 
                            cluster_columns = FALSE,
                            cluster_rows = FALSE, 
                            # cluster_rows = FALSE,
                            # show_row_names = FALSE,
                            column_names_side = 'top', 
                            rect_gp = gpar(col = 'white', lwd = 1), 
                            show_heatmap_legend = FALSE
                            # column_split = column_split_overlap,
                            # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D7high_rld_counts_2D_mCell_3D_mean,
                            #                                      gp = gpar(fontsize = 12))    )
)


# Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv <- melt(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv, id.vars = 'group')
# Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound <- melt(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound, id.vars = 'group')
# head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv)
# tail(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv)
# tail(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound)
# Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
#                                                            levels = c('0', '1', '2', '4'))
# Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
#                                                            levels = c('0', '1', '2', '4'))
# colnames(Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# colnames(Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# 
# ggplot(data = Gene_Cell_signature_D7high_rld_counts_2D_mCell_solv, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell recovery signature - 2D without compound') + 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_2D_solv_recovery_signature.tiff', 
#   dpi = 300
# )
# 
# ggplot(data = Gene_Cell_signature_D7high_rld_counts_2D_mCell_compound, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell recovery signature - 2D with compound')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_2D_compound_recovery_signature.tiff', 
#   dpi = 300
# )
# 
# 
# # invivo Gene signatures D7 high
# Gene_Cell_signature_D7high_rld_counts_invivo <- melt(Gene_Cell_signature_D7high_rld_counts_invivo, id.vars = 'group')
# Gene_Cell_signature_D7high_rld_counts_invivo$group <- factor(Gene_Cell_signature_D7high_rld_counts_invivo$group, 
#                                                     levels = c('0', '1', '3', '7'))
# colnames(Gene_Cell_signature_D7high_rld_counts_invivo) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# ggplot(data = Gene_Cell_signature_D7high_rld_counts_invivo, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell recovery signature - in vivo acute injury')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_invivo_acute_injury_recovery_signature.tiff', 
#   dpi = 300
# )
# 
# 
# # 3D Gene signatures D7 high
# Gene_Cell_signature_D7high_rld_counts_3D <- melt(Gene_Cell_signature_D7high_rld_counts_3D, id.vars = 'group')
# Gene_Cell_signature_D7high_rld_counts_3D$group <- factor(Gene_Cell_signature_D7high_rld_counts_3D$group, 
#                                                       levels = c('0', '4', '7', '10', '14'))
# colnames(Gene_Cell_signature_D7high_rld_counts_3D) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# ggplot(data = Gene_Cell_signature_D7high_rld_counts_3D, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell recovery signature - 3D')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_3D_recovery_signature.tiff', 
#   dpi = 300
# )




# Gene signatures D1 high - damage signature =======
Gene_Cell_signature_D1high_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                               Gene_Cell_signature_D1high, ]
# Gene_Cell_signature_D1high_rld_counts_2D_mCell_invitro <- as.data.frame(t(Gene_Cell_signature_D1high_rld_counts_2D_mCell_invitro))
# Gene_Cell_signature_D1high_rld_counts_2D_mCell_invitro$group <- c(rep('D0', 4), rep('D1_Solv', 4), 
#                                                                    rep('D2_Solv', 4), rep('D4_Solv',4), 
#                                                                    rep('D1_compound', 4), 
#                                                                    rep('D2_compound', 4), rep('D4_compound',4)
#                                                                   )
head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_invitro)

Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv <- Gene_Cell_signature_D1high_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv)

Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound <- Gene_Cell_signature_D1high_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound)


Gene_Cell_signature_D1high_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                    Gene_Cell_signature_D1high, ]
# Gene_Cell_signature_D1high_rld_counts_invivo <- as.data.frame(t(Gene_Cell_signature_D1high_rld_counts_invivo))
# Gene_Cell_signature_D1high_rld_counts_invivo$group <- c(rep('0', 4), rep('1', 2), rep('3', 4), rep('7',4))
head(Gene_Cell_signature_D1high_rld_counts_invivo)


Gene_Cell_signature_D1high_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                        Gene_Cell_signature_D1high, ]
# Gene_Cell_signature_D1high_rld_counts_3D <- as.data.frame(t(Gene_Cell_signature_D1high_rld_counts_3D))
# Gene_Cell_signature_D1high_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_Cell_signature_D1high_rld_counts_3D)


library(ggplot2)
library(reshape2)

# invitro Gene signatures D7 high

Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv))))
head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv)

Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv, 
                                                                   `D0` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean, 
                                                                   `D1` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean, 
                                                                   `D2` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean, 
                                                                   `D4` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean <- Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean)



Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound))))
head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound)

Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound, 
                                                                   `D0` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean, 
                                                                   `D1` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean, 
                                                                   `D2` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean, 
                                                                   `D4` = apply(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean <- Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean)



Gene_Cell_signature_D1high_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_Cell_signature_D1high_rld_counts_invivo))))
head(Gene_Cell_signature_D1high_rld_counts_invivo)

Gene_Cell_signature_D1high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_invivo, 
                                                            `D0` = apply(Gene_Cell_signature_D1high_rld_counts_invivo[1:4], 1, mean))
Gene_Cell_signature_D1high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_invivo_mean, 
                                                            `D1` = apply(Gene_Cell_signature_D1high_rld_counts_invivo[5: 6], 1, mean))
Gene_Cell_signature_D1high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_invivo_mean, 
                                                            `D3` = apply(Gene_Cell_signature_D1high_rld_counts_invivo[7: 10], 1, mean))
Gene_Cell_signature_D1high_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_invivo_mean, 
                                                            `D7` = apply(Gene_Cell_signature_D1high_rld_counts_invivo[11: 14], 1, mean))

Gene_Cell_signature_D1high_rld_counts_invivo_mean <- Gene_Cell_signature_D1high_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_Cell_signature_D1high_rld_counts_invivo_mean)


Gene_Cell_signature_D1high_rld_counts_3D <- as.data.frame(t(scale(t(Gene_Cell_signature_D1high_rld_counts_3D))))

Gene_Cell_signature_D1high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_3D, 
                                                              `D0` = apply(Gene_Cell_signature_D1high_rld_counts_3D[1: 3], 1, mean))
Gene_Cell_signature_D1high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_3D_mean, 
                                                              `D4` = apply(Gene_Cell_signature_D1high_rld_counts_3D_mean[4: 6], 1, mean))
Gene_Cell_signature_D1high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_3D_mean, 
                                                              `D7` = apply(Gene_Cell_signature_D1high_rld_counts_3D_mean[7: 9], 1, mean))
Gene_Cell_signature_D1high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_3D_mean, 
                                                              `D10` = apply(Gene_Cell_signature_D1high_rld_counts_3D_mean[10: 12], 1, mean))
Gene_Cell_signature_D1high_rld_counts_3D_mean <- mutate(Gene_Cell_signature_D1high_rld_counts_3D_mean, 
                                                              `D14` = apply(Gene_Cell_signature_D1high_rld_counts_3D_mean[13: 15], 1, mean))
Gene_Cell_signature_D1high_rld_counts_3D_mean <- Gene_Cell_signature_D1high_rld_counts_3D_mean[,16: 20]
head(Gene_Cell_signature_D1high_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_D1high_solv <- Heatmap(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean,
                          col = color_ht,
                          width = unit(8*ncol(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean), 'mm'), 
                          height = unit(6* nrow(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean), 'mm'),
                          row_title = 'Damaged gene signature', 
                          cluster_columns = FALSE,
                          cluster_rows = FALSE, 
                          # cluster_rows = FALSE,
                          # show_row_names = FALSE,
                          column_names_side = 'top', 
                          rect_gp = gpar(col = 'white', lwd = 1), 
                          show_heatmap_legend = FALSE
                          # column_split = column_split_overlap,
                          # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv_mean,
                          #                                      gp = gpar(fontsize = 12))    )
)

ht_D1high_compound <- Heatmap(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean,
                          col = color_ht,
                          width = unit(8*ncol(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean), 'mm'), 
                          height = unit(6* nrow(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean), 'mm'),
                          row_title = 'Damaged gene signature', 
                          cluster_columns = FALSE,
                          cluster_rows = FALSE, 
                          # cluster_rows = FALSE,
                          # show_row_names = FALSE,
                          column_names_side = 'top', 
                          rect_gp = gpar(col = 'white', lwd = 1), 
                          show_heatmap_legend = FALSE
                          # column_split = column_split_overlap,
                          # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound_mean,
                          #                                      gp = gpar(fontsize = 12))    )
)

ht_D1high_invivo <- Heatmap(Gene_Cell_signature_D1high_rld_counts_invivo_mean,
                            col = color_ht,
                            width = unit(8*ncol(Gene_Cell_signature_D1high_rld_counts_invivo_mean), 'mm'), 
                            height = unit(6* nrow(Gene_Cell_signature_D1high_rld_counts_invivo_mean), 'mm'),
                            row_title = 'Damaged gene signature', 
                            cluster_columns = FALSE,
                            cluster_rows = FALSE, 
                            # cluster_rows = FALSE,
                            # show_row_names = FALSE,
                            column_names_side = 'top', 
                            rect_gp = gpar(col = 'white', lwd = 1), 
                            show_heatmap_legend = FALSE
                            # column_split = column_split_overlap,
                            # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D1high_rld_counts_2D_mCell_invivo_mean,
                            #                                      gp = gpar(fontsize = 12))    )
)

ht_D1high_3D <- Heatmap(Gene_Cell_signature_D1high_rld_counts_3D_mean,
                              col = color_ht,
                              width = unit(8*ncol(Gene_Cell_signature_D1high_rld_counts_3D_mean), 'mm'), 
                              height = unit(6* nrow(Gene_Cell_signature_D1high_rld_counts_3D_mean), 'mm'),
                              row_title = 'Damaged gene signature', 
                              cluster_columns = FALSE,
                              cluster_rows = FALSE, 
                              # cluster_rows = FALSE,
                              # show_row_names = FALSE,
                              column_names_side = 'top', 
                              rect_gp = gpar(col = 'white', lwd = 1), 
                              show_heatmap_legend = FALSE
                              # column_split = column_split_overlap,
                              # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_D1high_rld_counts_2D_mCell_3D_mean,
                              #                                      gp = gpar(fontsize = 12))    )
)


# 
# library(ggplot2)
# library(reshape2)

# invitro Gene signatures D1 high
# Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv <- melt(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv, id.vars = 'group')
# Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound <- melt(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound, id.vars = 'group')
# head(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv)
# tail(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv)
# tail(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound)
# Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
#                                                                     levels = c('0', '1', '2', '4'))
# Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
#                                                                     levels = c('0', '1', '2', '4'))
# colnames(Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# colnames(Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# 
# ggplot(data = Gene_Cell_signature_D1high_rld_counts_2D_mCell_solv, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell damage signature - 2D without compound') + 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_2D_solv_damage_signature.tiff', 
#   dpi = 300
# )
# 
# ggplot(data = Gene_Cell_signature_D1high_rld_counts_2D_mCell_compound, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell damage signature - 2D with compound')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_2D_compound_damage_signature.tiff', 
#   dpi = 300
# )
# 
# 
# 
# # invivo Gene signatures D1 high
# Gene_Cell_signature_D1high_rld_counts_invivo <- melt(Gene_Cell_signature_D1high_rld_counts_invivo, id.vars = 'group')
# Gene_Cell_signature_D1high_rld_counts_invivo$group <- factor(Gene_Cell_signature_D1high_rld_counts_invivo$group, 
#                                                              levels = c('0', '1', '3', '7'))
# colnames(Gene_Cell_signature_D1high_rld_counts_invivo) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# ggplot(data = Gene_Cell_signature_D1high_rld_counts_invivo, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell damage signature - in vivo acute injury')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_invivo_damage_signature.tiff', 
#   dpi = 300
# )
# 
# 
# # 3D Gene signatures D1 high
# Gene_Cell_signature_D1high_rld_counts_3D <- melt(Gene_Cell_signature_D1high_rld_counts_3D, id.vars = 'group')
# Gene_Cell_signature_D1high_rld_counts_3D$group <- factor(Gene_Cell_signature_D1high_rld_counts_3D$group, 
#                                                                levels = c('0', '4', '7', '10', '14'))
# colnames(Gene_Cell_signature_D1high_rld_counts_3D) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# ggplot(data = Gene_Cell_signature_D1high_rld_counts_3D, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell damage signature - 3D')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_3D_damage_signature.tiff', 
#   dpi = 300
# )



# Gene signatures Healthy high - Healthy signature =====

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                               Gene_Cell_signature_Healthyhigh, ]
# Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invitro <- as.data.frame(t(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invitro))
# Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invitro$group <- c(rep('D0', 4), rep('D1_Solv', 4), 
#                                                                    rep('D2_Solv', 4), rep('D4_Solv',4), 
#                                                                    rep('D1_compound', 4), 
#                                                                    rep('D2_compound', 4), rep('D4_compound',4)
#                                                                   )
head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invitro)

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv <- Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv)

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound <- Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound)


Gene_Cell_signature_Healthyhigh_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                    Gene_Cell_signature_Healthyhigh, ]
# Gene_Cell_signature_Healthyhigh_rld_counts_invivo <- as.data.frame(t(Gene_Cell_signature_Healthyhigh_rld_counts_invivo))
# Gene_Cell_signature_Healthyhigh_rld_counts_invivo$group <- c(rep('0', 4), rep('1', 2), rep('3', 4), rep('7',4))
head(Gene_Cell_signature_Healthyhigh_rld_counts_invivo)


Gene_Cell_signature_Healthyhigh_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                        Gene_Cell_signature_Healthyhigh, ]
# Gene_Cell_signature_Healthyhigh_rld_counts_3D <- as.data.frame(t(Gene_Cell_signature_Healthyhigh_rld_counts_3D))
# Gene_Cell_signature_Healthyhigh_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_Cell_signature_Healthyhigh_rld_counts_3D)


library(ggplot2)
library(reshape2)

# invitro Gene signatures D7 high

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv))))
head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv)

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv, 
                                                                   `D0` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean, 
                                                                   `D1` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean, 
                                                                   `D2` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean, 
                                                                   `D4` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean <- Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean)



Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound))))
head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound)

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound, 
                                                                   `D0` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean, 
                                                                   `D1` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean, 
                                                                   `D2` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean, 
                                                                   `D4` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean <- Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean)



Gene_Cell_signature_Healthyhigh_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_Cell_signature_Healthyhigh_rld_counts_invivo))))
head(Gene_Cell_signature_Healthyhigh_rld_counts_invivo)

Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_invivo, 
                                                            `D0` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_invivo[1:4], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean, 
                                                            `D1` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_invivo[5: 6], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean, 
                                                            `D3` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_invivo[7: 10], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean, 
                                                            `D7` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_invivo[11: 14], 1, mean))

Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean <- Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean)


Gene_Cell_signature_Healthyhigh_rld_counts_3D <- as.data.frame(t(scale(t(Gene_Cell_signature_Healthyhigh_rld_counts_3D))))

Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_3D, 
                                                              `D0` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_3D[1: 3], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean, 
                                                              `D4` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean[4: 6], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean, 
                                                              `D7` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean[7: 9], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean, 
                                                              `D10` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean[10: 12], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean <- mutate(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean, 
                                                              `D14` = apply(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean[13: 15], 1, mean))
Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean <- Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean[,16: 20]
head(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_Healthyhigh_solv <- Heatmap(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean,
                          col = color_ht,
                          width = unit(8*ncol(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean), 'mm'), 
                          height = unit(6* nrow(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean), 'mm'),
                          row_title = 'Healthy gene signature', 
                          cluster_columns = FALSE,
                          cluster_rows = FALSE, 
                          # cluster_rows = FALSE,
                          # show_row_names = FALSE,
                          column_names_side = 'top', 
                          rect_gp = gpar(col = 'white', lwd = 1), 
                          heatmap_legend_param = list(title = 'Normalized counts')
                          # column_split = column_split_overlap,
                          # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv_mean,
                          #                                      gp = gpar(fontsize = 12))    )
)

ht_Healthyhigh_compound <- Heatmap(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean,
                          col = color_ht,
                          width = unit(8*ncol(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean), 'mm'), 
                          height = unit(6* nrow(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean), 'mm'),
                          row_title = 'Healthy gene signature', 
                          cluster_columns = FALSE,
                          cluster_rows = FALSE, 
                          # cluster_rows = FALSE,
                          # show_row_names = FALSE,
                          column_names_side = 'top', 
                          rect_gp = gpar(col = 'white', lwd = 1), 
                          heatmap_legend_param = list(title = 'Normalized counts')
                          # column_split = column_split_overlap,
                          # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound_mean,
                          #                                      gp = gpar(fontsize = 12))    )
)

ht_Healthyhigh_invivo <- Heatmap(Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean,
                            col = color_ht,
                            width = unit(8*ncol(Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean), 'mm'), 
                            height = unit(6* nrow(Gene_Cell_signature_Healthyhigh_rld_counts_invivo_mean), 'mm'),
                            row_title = 'Healthy gene signature', 
                            cluster_columns = FALSE,
                            cluster_rows = FALSE, 
                            # cluster_rows = FALSE,
                            # show_row_names = FALSE,
                            column_names_side = 'top', 
                            rect_gp = gpar(col = 'white', lwd = 1), 
                            heatmap_legend_param = list(title = 'Normalized counts')
                            # column_split = column_split_overlap,
                            # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_invivo_mean,
                            #                                      gp = gpar(fontsize = 12))    )
)

ht_Healthyhigh_3D <- Heatmap(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean,
                              col = color_ht,
                              width = unit(8*ncol(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean), 'mm'), 
                              height = unit(6* nrow(Gene_Cell_signature_Healthyhigh_rld_counts_3D_mean), 'mm'),
                              row_title = 'Healthy gene signature', 
                              cluster_columns = FALSE,
                              cluster_rows = FALSE, 
                              # cluster_rows = FALSE,
                              # show_row_names = FALSE,
                              column_names_side = 'top', 
                              rect_gp = gpar(col = 'white', lwd = 1), 
                              heatmap_legend_param = list(title = 'Normalized counts')
                              # column_split = column_split_overlap,
                              # row_names_max_width = max_text_width(rownames(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_3D_mean,
                              #                                      gp = gpar(fontsize = 12))    )
)


# library(ggplot2)
# library(reshape2)
# 
# # invitro Gene signatures Healthy high
# Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv <- melt(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv, id.vars = 'group')
# Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound <- melt(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound, id.vars = 'group')
# head(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv)
# tail(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv)
# tail(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound)
# Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
#                                                                     levels = c('0', '1', '2', '4'))
# Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound$group <- factor(c(rep('0', 4), rep('1',4), rep('2', 4), rep('4', 4)), 
#                                                                     levels = c('0', '1', '2', '4'))
# colnames(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# colnames(Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# 
# ggplot(data = Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_solv, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell healthy signature - 2D without compound') + 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_2D_solv_healthy_signature.tiff', 
#   dpi = 300
# )
# 
# ggplot(data = Gene_Cell_signature_Healthyhigh_rld_counts_2D_mCell_compound, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) + 
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell healthy signature - 2D with compound')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_2D_compound_healthy_signature.tiff', 
#   dpi = 300
# )
# 
# 
# 
# # invivo Gene signatures Healthy high
# Gene_Cell_signature_Healthyhigh_rld_counts_invivo <- melt(Gene_Cell_signature_Healthyhigh_rld_counts_invivo, id.vars = 'group')
# Gene_Cell_signature_Healthyhigh_rld_counts_invivo$group <- factor(Gene_Cell_signature_Healthyhigh_rld_counts_invivo$group, 
#                                                              levels = c('0', '1', '3', '7'))
# colnames(Gene_Cell_signature_Healthyhigh_rld_counts_invivo) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# ggplot(data = Gene_Cell_signature_Healthyhigh_rld_counts_invivo, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell Healthy signature - in vivo acute injury')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_Cell_invivo_healthy_signature.tiff', 
#   dpi = 300
# )
# 
# 
# # 3D Gene signatures Healthy high
# Gene_Cell_signature_Healthyhigh_rld_counts_3D <- melt(Gene_Cell_signature_Healthyhigh_rld_counts_3D, id.vars = 'group')
# Gene_Cell_signature_Healthyhigh_rld_counts_3D$group <- factor(Gene_Cell_signature_Healthyhigh_rld_counts_3D$group, 
#                                                                levels = c('0', '4', '7', '10', '14'))
# colnames(Gene_Cell_signature_Healthyhigh_rld_counts_3D) <- c('Timepoint', 'variable', 'log2(normalized_counts)')
# ggplot(data = Gene_Cell_signature_Healthyhigh_rld_counts_3D, aes(x = `Timepoint`, y = `log2(normalized_counts)`)) +
#   geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + ggtitle('Cell Healthy signature - 3D')+ 
#   geom_jitter(shape = 16, position = position_jitter(0.2))
# 
# setwd('D:/foldername/Cell_Dedif/20230906/plots/')
# ggplot2::ggsave(
#   filename = 'boxplot_3D_healthy_signature.tiff', 
#   dpi = 300
# )

# Gene signatures 4 damaged signature -----------
Gene_Cell_damage_signature_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                                    Gene_Cell_damage_signature, ]
# Gene_Cell_damage_signature_rld_counts_2D_mCell_invitro <- as.data.frame(t(Gene_Cell_damage_signature_rld_counts_2D_mCell_invitro))
# Gene_Cell_damage_signature_rld_counts_2D_mCell_invitro$group <- c(rep('D0', 4), rep('D1_Solv', 4), 
#                                                                    rep('D2_Solv', 4), rep('D4_Solv',4), 
#                                                                    rep('D1_compound', 4), 
#                                                                    rep('D2_compound', 4), rep('D4_compound',4)
#                                                                   )
head(Gene_Cell_damage_signature_rld_counts_2D_mCell_invitro)

Gene_Cell_damage_signature_rld_counts_2D_mCell_solv <- Gene_Cell_damage_signature_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv)

Gene_Cell_damage_signature_rld_counts_2D_mCell_compound <- Gene_Cell_damage_signature_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound)


Gene_Cell_damage_signature_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                         Gene_Cell_damage_signature, ]
# Gene_Cell_damage_signature_rld_counts_invivo <- as.data.frame(t(Gene_Cell_damage_signature_rld_counts_invivo))
# Gene_Cell_damage_signature_rld_counts_invivo$group <- c(rep('0', 4), rep('1', 2), rep('3', 4), rep('7',4))
head(Gene_Cell_damage_signature_rld_counts_invivo)


Gene_Cell_damage_signature_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                             Gene_Cell_damage_signature, ]
# Gene_Cell_damage_signature_rld_counts_3D <- as.data.frame(t(Gene_Cell_damage_signature_rld_counts_3D))
# Gene_Cell_damage_signature_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_Cell_damage_signature_rld_counts_3D)


library(ggplot2)
library(reshape2)


Gene_Cell_damage_signature_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv))))
head(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv)

Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv, 
                                                                        `D0` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean, 
                                                                        `D1` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean, 
                                                                        `D2` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean, 
                                                                        `D4` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean <- Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean)



Gene_Cell_damage_signature_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound))))
head(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound)

Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound, 
                                                                        `D0` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean, 
                                                                        `D1` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean, 
                                                                        `D2` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean, 
                                                                        `D4` = apply(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean <- Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean)



Gene_Cell_damage_signature_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_Cell_damage_signature_rld_counts_invivo))))
head(Gene_Cell_damage_signature_rld_counts_invivo)

Gene_Cell_damage_signature_rld_counts_invivo_mean <- mutate(Gene_Cell_damage_signature_rld_counts_invivo, 
                                                                 `D0` = apply(Gene_Cell_damage_signature_rld_counts_invivo[1:4], 1, mean))
Gene_Cell_damage_signature_rld_counts_invivo_mean <- mutate(Gene_Cell_damage_signature_rld_counts_invivo_mean, 
                                                                 `D1` = apply(Gene_Cell_damage_signature_rld_counts_invivo[5: 6], 1, mean))
Gene_Cell_damage_signature_rld_counts_invivo_mean <- mutate(Gene_Cell_damage_signature_rld_counts_invivo_mean, 
                                                                 `D3` = apply(Gene_Cell_damage_signature_rld_counts_invivo[7: 10], 1, mean))
Gene_Cell_damage_signature_rld_counts_invivo_mean <- mutate(Gene_Cell_damage_signature_rld_counts_invivo_mean, 
                                                                 `D7` = apply(Gene_Cell_damage_signature_rld_counts_invivo[11: 14], 1, mean))

Gene_Cell_damage_signature_rld_counts_invivo_mean <- Gene_Cell_damage_signature_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_Cell_damage_signature_rld_counts_invivo_mean)


Gene_Cell_damage_signature_rld_counts_3D <- as.data.frame(t(scale(t(Gene_Cell_damage_signature_rld_counts_3D))))

Gene_Cell_damage_signature_rld_counts_3D_mean <- mutate(Gene_Cell_damage_signature_rld_counts_3D, 
                                                                   `D0` = apply(Gene_Cell_damage_signature_rld_counts_3D[1: 3], 1, mean))
Gene_Cell_damage_signature_rld_counts_3D_mean <- mutate(Gene_Cell_damage_signature_rld_counts_3D_mean, 
                                                                   `D4` = apply(Gene_Cell_damage_signature_rld_counts_3D_mean[4: 6], 1, mean))
Gene_Cell_damage_signature_rld_counts_3D_mean <- mutate(Gene_Cell_damage_signature_rld_counts_3D_mean, 
                                                                   `D7` = apply(Gene_Cell_damage_signature_rld_counts_3D_mean[7: 9], 1, mean))
Gene_Cell_damage_signature_rld_counts_3D_mean <- mutate(Gene_Cell_damage_signature_rld_counts_3D_mean, 
                                                                   `D10` = apply(Gene_Cell_damage_signature_rld_counts_3D_mean[10: 12], 1, mean))
Gene_Cell_damage_signature_rld_counts_3D_mean <- mutate(Gene_Cell_damage_signature_rld_counts_3D_mean, 
                                                                   `D14` = apply(Gene_Cell_damage_signature_rld_counts_3D_mean[13: 15], 1, mean))
Gene_Cell_damage_signature_rld_counts_3D_mean <- Gene_Cell_damage_signature_rld_counts_3D_mean[,16: 20]
head(Gene_Cell_damage_signature_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_damage_signature_solv <- Heatmap(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean,
                               col = color_ht,
                               width = unit(8*ncol(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean), 'mm'), 
                               height = unit(6* nrow(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean), 'mm'),
                               column_title = 'In vitro 2D without compound',
                               cluster_columns = FALSE,
                               cluster_rows = FALSE, 
                               # cluster_rows = FALSE,
                               # show_row_names = FALSE,
                               column_names_side = 'top', 
                               rect_gp = gpar(col = 'white', lwd = 1), 
                               heatmap_legend_param = list(title = 'Normalized counts')
                               # column_split = column_split_overlap,
                               # row_names_max_width = max_text_width(rownames(Gene_Cell_damage_signature_rld_counts_2D_mCell_solv_mean,
                               #                                      gp = gpar(fontsize = 12))    )
)

ht_damage_signature_compound <- Heatmap(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean,
                               col = color_ht,
                               width = unit(8*ncol(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean), 'mm'), 
                               height = unit(6* nrow(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean), 'mm'),
                               column_title = 'In vitro 2D with compound',
                               cluster_columns = FALSE,
                               cluster_rows = FALSE, 
                               # cluster_rows = FALSE,
                               # show_row_names = FALSE,
                               column_names_side = 'top', 
                               rect_gp = gpar(col = 'white', lwd = 1), 
                               heatmap_legend_param = list(title = 'Normalized counts')
                               # column_split = column_split_overlap,
                               # row_names_max_width = max_text_width(rownames(Gene_Cell_damage_signature_rld_counts_2D_mCell_compound_mean,
                               #                                      gp = gpar(fontsize = 12))    )
)

ht_damage_signature_invivo <- Heatmap(Gene_Cell_damage_signature_rld_counts_invivo_mean,
                                 col = color_ht,
                                 width = unit(8*ncol(Gene_Cell_damage_signature_rld_counts_invivo_mean), 'mm'), 
                                 height = unit(6* nrow(Gene_Cell_damage_signature_rld_counts_invivo_mean), 'mm'),
                                 column_title = 'In vivo acute injury',
                                 cluster_columns = FALSE,
                                 cluster_rows = FALSE, 
                                 # cluster_rows = FALSE,
                                 # show_row_names = FALSE,
                                 column_names_side = 'top', 
                                 rect_gp = gpar(col = 'white', lwd = 1), 
                                 heatmap_legend_param = list(title = 'Normalized counts')
                                 # column_split = column_split_overlap,
                                 # row_names_max_width = max_text_width(rownames(Gene_Cell_damage_signature_rld_counts_2D_mCell_invivo_mean,
                                 #                                      gp = gpar(fontsize = 12))    )
)

ht_damage_signature_3D <- Heatmap(Gene_Cell_damage_signature_rld_counts_3D_mean,
                                   col = color_ht,
                                   width = unit(8*ncol(Gene_Cell_damage_signature_rld_counts_3D_mean), 'mm'), 
                                   height = unit(6* nrow(Gene_Cell_damage_signature_rld_counts_3D_mean), 'mm'),
                                   column_title = 'Moues 3D',
                                   cluster_columns = FALSE,
                                   cluster_rows = FALSE, 
                                   # cluster_rows = FALSE,
                                   # show_row_names = FALSE,
                                   column_names_side = 'top', 
                                   rect_gp = gpar(col = 'white', lwd = 1), 
                                   heatmap_legend_param = list(title = 'Normalized counts')
                                   # column_split = column_split_overlap,
                                   # row_names_max_width = max_text_width(rownames(Gene_Cell_damage_signature_rld_counts_2D_mCell_3D_mean,
                                   #                                      gp = gpar(fontsize = 12))    )
)
# Heatmap of gene signatures-----
draw(ht_Healthyhigh_solv %v% ht_D1high_solv %v% ht_D7high_solv)
draw(ht_Healthyhigh_compound %v% ht_D1high_compound %v% ht_D7high_compound)
draw(ht_Healthyhigh_invivo %v% ht_D1high_invivo %v% ht_D7high_invivo)
draw(ht_Healthyhigh_3D %v% ht_D1high_3D %v% ht_D7high_3D)
# draw it separately
# draw(ht_damage_signature_solv + ht_damage_signature_compound + ht_damage_signature_invivo + ht_damage_signature_3D)
draw(ht_damage_signature_solv)
draw(ht_damage_signature_compound)
draw(ht_damage_signature_invivo)
draw(ht_damage_signature_3D)

# DE genes between solv and compound -----
test1 <- DE_2D_mCell_timepoint_compound_D1[!(rownames(DE_2D_mCell_timepoint_compound_D1) %in% rownames(DE_2D_mCell_timepoint_solv_D1)), ]
test2_compound <- DE_2D_mCell_timepoint_compound_D2[!(rownames(DE_2D_mCell_timepoint_compound_D2) %in% rownames(DE_2D_mCell_timepoint_solv_D2)), ]
test2_solv <- DE_2D_mCell_timepoint_solv_D2[!(rownames(DE_2D_mCell_timepoint_solv_D2) %in% rownames(DE_2D_mCell_timepoint_compound_D2)), ]
test3_compound <- DE_2D_mCell_timepoint_compound_D4[!(rownames(DE_2D_mCell_timepoint_compound_D4) %in% rownames(DE_2D_mCell_timepoint_solv_D4)), ]
test3_solv <- DE_2D_mCell_timepoint_solv_D4[!(rownames(DE_2D_mCell_timepoint_solv_D4) %in% rownames(DE_2D_mCell_timepoint_compound_D4)), ]




# In the chronic injury(2w, 4w from scRNA data) the DE gene expression in the 2D culture------------
DEgene_chronic_injury_4w <- c('Fabp4', 'Rps18-ps6', 'Egr1', 'Cd74', 'Hp', 'Ednrb', 'H2-Eb1', 'H2-Aa')
Gene_chronic_injury_4w_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                               DEgene_chronic_injury_4w, ]
# Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro <- as.data.frame(t(Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro))
# Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro$group <- c(rep('D0', 4), rep('D1_Solv', 4), 
#                                                                    rep('D2_Solv', 4), rep('D4_Solv',4), 
#                                                                    rep('D1_compound', 4), 
#                                                                    rep('D2_compound', 4), rep('D4_compound',4)
#                                                                   )
head(Gene_chronic_injury_4w_rld_counts_2D_mCell_invitro)

Gene_chronic_injury_4w_rld_counts_2D_mCell_solv <- Gene_chronic_injury_4w_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv)

Gene_chronic_injury_4w_rld_counts_2D_mCell_compound <- Gene_chronic_injury_4w_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound)


Gene_chronic_injury_4w_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                DEgene_chronic_injury_4w, ]
# Gene_chronic_injury_4w_rld_counts_invivo <- as.data.frame(t(Gene_chronic_injury_4w_rld_counts_invivo))
# Gene_chronic_injury_4w_rld_counts_invivo$group <- c(rep('0', 4), rep('1', 2), rep('3', 4), rep('7',4))
head(Gene_chronic_injury_4w_rld_counts_invivo)


Gene_chronic_injury_4w_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                    DEgene_chronic_injury_4w, ]
# Gene_chronic_injury_4w_rld_counts_3D <- as.data.frame(t(Gene_chronic_injury_4w_rld_counts_3D))
# Gene_chronic_injury_4w_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_chronic_injury_4w_rld_counts_3D)


Gene_chronic_injury_4w_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv))))
head(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv)

Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv, 
                                                                   `D0` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean, 
                                                                   `D1` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean, 
                                                                   `D2` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean, 
                                                                   `D4` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean <- Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean)



Gene_chronic_injury_4w_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound))))
head(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound)

Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound, 
                                                                   `D0` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean, 
                                                                   `D1` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean, 
                                                                   `D2` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean <- mutate(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean, 
                                                                   `D4` = apply(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean <- Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean)



Gene_chronic_injury_4w_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_chronic_injury_4w_rld_counts_invivo))))
head(Gene_chronic_injury_4w_rld_counts_invivo)

Gene_chronic_injury_4w_rld_counts_invivo_mean <- mutate(Gene_chronic_injury_4w_rld_counts_invivo, 
                                                            `D0` = apply(Gene_chronic_injury_4w_rld_counts_invivo[1:4], 1, mean))
Gene_chronic_injury_4w_rld_counts_invivo_mean <- mutate(Gene_chronic_injury_4w_rld_counts_invivo_mean, 
                                                            `D1` = apply(Gene_chronic_injury_4w_rld_counts_invivo[5: 6], 1, mean))
Gene_chronic_injury_4w_rld_counts_invivo_mean <- mutate(Gene_chronic_injury_4w_rld_counts_invivo_mean, 
                                                            `D3` = apply(Gene_chronic_injury_4w_rld_counts_invivo[7: 10], 1, mean))
Gene_chronic_injury_4w_rld_counts_invivo_mean <- mutate(Gene_chronic_injury_4w_rld_counts_invivo_mean, 
                                                            `D7` = apply(Gene_chronic_injury_4w_rld_counts_invivo[11: 14], 1, mean))

Gene_chronic_injury_4w_rld_counts_invivo_mean <- Gene_chronic_injury_4w_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_chronic_injury_4w_rld_counts_invivo_mean)


Gene_chronic_injury_4w_rld_counts_3D <- as.data.frame(t(scale(t(Gene_chronic_injury_4w_rld_counts_3D))))

Gene_chronic_injury_4w_rld_counts_3D_mean <- mutate(Gene_chronic_injury_4w_rld_counts_3D, 
                                                              `D0` = apply(Gene_chronic_injury_4w_rld_counts_3D[1: 3], 1, mean))
Gene_chronic_injury_4w_rld_counts_3D_mean <- mutate(Gene_chronic_injury_4w_rld_counts_3D_mean, 
                                                              `D4` = apply(Gene_chronic_injury_4w_rld_counts_3D_mean[4: 6], 1, mean))
Gene_chronic_injury_4w_rld_counts_3D_mean <- mutate(Gene_chronic_injury_4w_rld_counts_3D_mean, 
                                                              `D7` = apply(Gene_chronic_injury_4w_rld_counts_3D_mean[7: 9], 1, mean))
Gene_chronic_injury_4w_rld_counts_3D_mean <- mutate(Gene_chronic_injury_4w_rld_counts_3D_mean, 
                                                              `D10` = apply(Gene_chronic_injury_4w_rld_counts_3D_mean[10: 12], 1, mean))
Gene_chronic_injury_4w_rld_counts_3D_mean <- mutate(Gene_chronic_injury_4w_rld_counts_3D_mean, 
                                                              `D14` = apply(Gene_chronic_injury_4w_rld_counts_3D_mean[13: 15], 1, mean))
Gene_chronic_injury_4w_rld_counts_3D_mean <- Gene_chronic_injury_4w_rld_counts_3D_mean[,16: 20]
head(Gene_chronic_injury_4w_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_chronic_injury_4w_solv <- Heatmap(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean,
                                    col = color_ht,
                                    width = unit(8*ncol(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean), 'mm'), 
                                    height = unit(6* nrow(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean), 'mm'),
                                    column_title = 'In vitro 2D without compound',
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE, 
                                    # cluster_rows = FALSE,
                                    # show_row_names = FALSE,
                                    column_names_side = 'top', 
                                    rect_gp = gpar(col = 'white', lwd = 1), 
                                    heatmap_legend_param = list(title = 'Normalized counts')
                                    # column_split = column_split_overlap,
                                    # row_names_max_width = max_text_width(rownames(Gene_chronic_injury_4w_rld_counts_2D_mCell_solv_mean,
                                    #                                      gp = gpar(fontsize = 12))    )
)

ht_chronic_injury_4w_compound <- Heatmap(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean,
                                    col = color_ht,
                                    width = unit(8*ncol(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean), 'mm'), 
                                    height = unit(6* nrow(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean), 'mm'),
                                    column_title = 'In vitro 2D with compound',
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE, 
                                    # cluster_rows = FALSE,
                                    # show_row_names = FALSE,
                                    column_names_side = 'top', 
                                    rect_gp = gpar(col = 'white', lwd = 1), 
                                    heatmap_legend_param = list(title = 'Normalized counts')
                                    # column_split = column_split_overlap,
                                    # row_names_max_width = max_text_width(rownames(Gene_chronic_injury_4w_rld_counts_2D_mCell_compound_mean,
                                    #                                      gp = gpar(fontsize = 12))    )
)

ht_chronic_injury_4w_invivo <- Heatmap(Gene_chronic_injury_4w_rld_counts_invivo_mean,
                                      col = color_ht,
                                      width = unit(8*ncol(Gene_chronic_injury_4w_rld_counts_invivo_mean), 'mm'), 
                                      height = unit(6* nrow(Gene_chronic_injury_4w_rld_counts_invivo_mean), 'mm'),
                                      column_title = 'In vivo acute injury',
                                      cluster_columns = FALSE,
                                      cluster_rows = FALSE, 
                                      # cluster_rows = FALSE,
                                      # show_row_names = FALSE,
                                      column_names_side = 'top', 
                                      rect_gp = gpar(col = 'white', lwd = 1), 
                                      heatmap_legend_param = list(title = 'Normalized counts')
                                      # column_split = column_split_overlap,
                                      # row_names_max_width = max_text_width(rownames(Gene_chronic_injury_4w_rld_counts_2D_mCell_invivo_mean,
                                      #                                      gp = gpar(fontsize = 12))    )
)

ht_chronic_injury_4w_3D <- Heatmap(Gene_chronic_injury_4w_rld_counts_3D_mean,
                                        col = color_ht,
                                        width = unit(8*ncol(Gene_chronic_injury_4w_rld_counts_3D_mean), 'mm'), 
                                        height = unit(6* nrow(Gene_chronic_injury_4w_rld_counts_3D_mean), 'mm'),
                                        column_title = 'Moues 3D',
                                        cluster_columns = FALSE,
                                        cluster_rows = FALSE, 
                                        # cluster_rows = FALSE,
                                        # show_row_names = FALSE,
                                        column_names_side = 'top', 
                                        rect_gp = gpar(col = 'white', lwd = 1), 
                                        heatmap_legend_param = list(title = 'Normalized counts')
                                        # column_split = column_split_overlap,
                                        # row_names_max_width = max_text_width(rownames(Gene_chronic_injury_4w_rld_counts_2D_mCell_3D_mean,
                                        #                                      gp = gpar(fontsize = 12))    )
)

draw(ht_chronic_injury_4w_solv)
draw(ht_chronic_injury_4w_compound)
draw(ht_chronic_injury_4w_invivo)
draw(ht_chronic_injury_4w_3D)





# In vitro, invivo and 3Ds signatures change - GSVA ========

# Cell gene signature GSVA===========
library(GSVA)
head(Gene_Cell_signature_D7high_rld_counts_2D_mCell_invitro)
Gene_Cell_signature_Recovery_list <- list("Recovery Cell signature" = Gene_Cell_signature_D7high)
head(Gene_Cell_signature_Recovery_list)
Gene_Cell_signature_Healthy_list <- list("Healthy Cell signature" = Gene_Cell_signature_Healthyhigh)
head(Gene_Cell_signature_Healthy_list)
Gene_Cell_signature_Damaged_list <- list("Damaged Cell signature" = c('Ccl6', 'Cd14', 'Ccl9', 'Insr', 'Adm', 'Nid1', 
                                                                      'Timp1', 'Akap12', 'Lyve1', 'Clec1b', 'F2r', 
                                                                      'Tfpi2', 'Cd4', 'Il6st', 'Crhbp', 'Cd36', 
                                                                      'Stab1','Fabp5','Vwf','Vwa1','Fabp4'))
head(Gene_Cell_signature_Damaged_list)
# chronic injury
Gene_chronic_injury_4w_list <- list("Chronic injury 4w signature" = DEgene_chronic_injury_4w)
head(Gene_chronic_injury_4w_list)

# Fisrt gsva
gsva_Recovery_invitro <- gsva(normcounts_2D_mCell_timepoint, Gene_Cell_signature_Recovery_list)
head(gsva_Recovery_invitro)
gsva_Recovery_solv <- as.matrix(as.data.frame(gsva_Recovery_invitro)[1:16])
head(gsva_Recovery_solv)
gsva_Recovery_compound <- as.matrix(as.data.frame(gsva_Recovery_invitro)[c(1:4, 17:28)])
head(gsva_Recovery_compound)
gsva_Recovery_invivo <- gsva(normcounts_invivo_mCell_timepoint_acute, Gene_Cell_signature_Recovery_list)
head(gsva_Recovery_invivo)
gsva_Recovery_3D <- gsva(normcounts_3D, Gene_Cell_signature_Recovery_list)
head(gsva_Recovery_3D)


gsva_Healthy_invitro <- gsva(normcounts_2D_mCell_timepoint, Gene_Cell_signature_Healthy_list)
head(gsva_Healthy_invitro)
gsva_Healthy_solv <- as.matrix(as.data.frame(gsva_Healthy_invitro)[1:16])
head(gsva_Healthy_solv)
gsva_Healthy_compound <- as.matrix(as.data.frame(gsva_Healthy_invitro)[c(1:4, 17:28)])
head(gsva_Healthy_compound)
gsva_Healthy_invivo <- gsva(normcounts_invivo_mCell_timepoint_acute, Gene_Cell_signature_Healthy_list)
head(gsva_Healthy_invivo)
gsva_Healthy_3D <- gsva(normcounts_3D, Gene_Cell_signature_Healthy_list)
head(gsva_Healthy_3D)

gsva_Damaged_invitro <- gsva(normcounts_2D_mCell_timepoint, Gene_Cell_signature_Damaged_list)
head(gsva_Damaged_invitro)
gsva_Damaged_solv <- as.matrix(as.data.frame(gsva_Damaged_invitro)[1:16])
head(gsva_Damaged_solv)
gsva_Damaged_compound <- as.matrix(as.data.frame(gsva_Damaged_invitro)[c(1:4, 17:28)])
head(gsva_Damaged_compound)
gsva_Damaged_invivo <- gsva(normcounts_invivo_mCell_timepoint_acute, Gene_Cell_signature_Damaged_list)
head(gsva_Damaged_invivo)
gsva_Damaged_3D <- gsva(normcounts_3D, Gene_Cell_signature_Damaged_list)
head(gsva_Damaged_3D)



gsva_solv <- rbind(gsva_Healthy_solv, gsva_Damaged_solv, gsva_Recovery_solv)
head(gsva_solv)

gsva_compound <- rbind(gsva_Healthy_compound, gsva_Damaged_compound, gsva_Recovery_compound)
head(gsva_compound)

gsva_invivo <- rbind(gsva_Healthy_invivo, gsva_Damaged_invivo, gsva_Recovery_invivo)
head(gsva_invivo)

gsva_3D <- rbind(gsva_Healthy_3D, gsva_Damaged_3D, gsva_Recovery_3D)
head(gsva_3D)


# calculate FDR by limma=====
library(limma)
# solv
mod_solv <- model.matrix(~ factor(c(rep('Day0', 4), rep('Day1', 4), rep('Day2', 4), rep('Day4', 4))))
colnames(mod_solv) <- c('Day0', 'Day1', 'Day2', 'Day4')
fit_solv <- lmFit(gsva_solv, mod_solv)
head(coef(fit_solv))
# Day1 - Day0
contr_solv <- makeContrasts(contrast = Day0-Day1,  levels = coef(fit_solv))
contr_solv
fit_solv <- contrasts.fit(fit_solv, contr_solv)
fit_solv <- eBayes(fit_solv)
tt_solv <- topTable(fit_solv)
tt_solv
FDR_gsva_solv_D1 <- tt_solv[, c('P.Value',"adj.P.Val")]
FDR_gsva_solv_D1$rn <- rownames(FDR_gsva_solv_D1)
colnames(FDR_gsva_solv_D1) <- c('P.Value.D1',"adj.P.Val.D1", 'rn')
head(FDR_gsva_solv_D1)
# Day2 - Day0
contr_solv <- makeContrasts(contrast = Day0-Day2, levels = colnames(coef(fit_solv)))
contr_solv
fit_solv <- contrasts.fit(fit_solv, contr_solv)
fit_solv <- eBayes(fit_solv)
tt_solv <- topTable(fit_solv)
tt_solv
FDR_gsva_solv_D2 <- tt_solv[, c('P.Value',"adj.P.Val")]
FDR_gsva_solv_D2$rn <- rownames(FDR_gsva_solv_D2)
colnames(FDR_gsva_solv_D2) <- c('P.Value.D2',"adj.P.Val.D2", 'rn')
# Day4 - Day0
contr_solv <- makeContrasts(contrast = Day0-Day4, levels = colnames(coef(fit_solv)))
contr_solv
fit_solv <- contrasts.fit(fit_solv, contr_solv)
fit_solv <- eBayes(fit_solv)
tt_solv <- topTable(fit_solv)
tt_solv
FDR_gsva_solv_D4 <- tt_solv[,c('P.Value',"adj.P.Val")]
FDR_gsva_solv_D4$rn <- rownames(FDR_gsva_solv_D4)
colnames(FDR_gsva_solv_D4) <- c('P.Value.D4',"adj.P.Val.D4", 'rn')
# combine FDRs
Multimerge <- function(x, y){
              df <- merge(x, y, by= "rn", all.x= TRUE, all.y= TRUE)
              return(df)
}
FDR_gsva_solv <- Reduce(Multimerge, list(FDR_gsva_solv_D1, FDR_gsva_solv_D2, FDR_gsva_solv_D4))
head(FDR_gsva_solv)
rownames(FDR_gsva_solv) <- FDR_gsva_solv$rn
FDR_gsva_solv <- FDR_gsva_solv[,-c(1, 2, 4, 6)]
head(FDR_gsva_solv)

#compound
mod_compound <- model.matrix(~ factor(c(rep('Day0', 4), rep('Day1', 4), rep('Day2', 4), rep('Day4', 4))))
colnames(mod_compound) <- c('Day0', 'Day1', 'Day2', 'Day4')
fit_compound <- lmFit(gsva_compound, mod_compound)
head(coef(fit_compound))
# Day1 - Day0
contr_compound <- makeContrasts(contrast = Day0-Day1, levels = colnames(coef(fit_compound)))
fit_compound <- contrasts.fit(fit_compound, contr_compound)
fit_compound <- eBayes(fit_compound)
tt_compound <- topTable(fit_compound)
tt_compound
FDR_gsva_compound_D1 <- tt_compound[, c('P.Value',"adj.P.Val")]
FDR_gsva_compound_D1$rn <- rownames(FDR_gsva_compound_D1)
colnames(FDR_gsva_compound_D1) <- c('P.Value.D1',"adj.P.Val.D1", 'rn')
head(FDR_gsva_compound_D1)
# Day2 - Day0
contr_compound <- makeContrasts(contrast = Day0-Day2, levels = colnames(coef(fit_compound)))
fit_compound <- contrasts.fit(fit_compound, contr_compound)
fit_compound <- eBayes(fit_compound)
tt_compound <- topTable(fit_compound)
tt_compound
FDR_gsva_compound_D2 <- tt_compound[, c('P.Value',"adj.P.Val")]
FDR_gsva_compound_D2$rn <- rownames(FDR_gsva_compound_D2)
colnames(FDR_gsva_compound_D2) <- c('P.Value.D2',"adj.P.Val.D2", 'rn')
# Day4 - Day0
contr_compound <- makeContrasts(contrast = Day0-Day4, levels = colnames(coef(fit_compound)))
fit_compound <- contrasts.fit(fit_compound, contr_compound)
fit_compound <- eBayes(fit_compound)
tt_compound <- topTable(fit_compound)
tt_compound
FDR_gsva_compound_D4 <- tt_compound[, c('P.Value',"adj.P.Val")]
FDR_gsva_compound_D4$rn <- rownames(FDR_gsva_compound_D4)
colnames(FDR_gsva_compound_D4) <- c('P.Value.D4',"adj.P.Val.D4", 'rn')
# combine FDRs
FDR_gsva_compound <- Reduce(Multimerge, list(FDR_gsva_compound_D1, FDR_gsva_compound_D2, FDR_gsva_compound_D4))
head(FDR_gsva_compound)
rownames(FDR_gsva_compound) <- FDR_gsva_compound$rn
FDR_gsva_compound <- FDR_gsva_compound[,-c(1, 2, 4, 6)]
head(FDR_gsva_compound)

#invivo
mod_invivo <- model.matrix(~ factor(c(rep('Day0', 4), rep('Day1', 2), rep('Day3', 4), rep('Day7', 4))))
colnames(mod_invivo) <- c('Day0', 'Day1', 'Day3', 'Day7')
fit_invivo <- lmFit(gsva_invivo, mod_invivo)
head(coef(fit_invivo))
# Day1 - Day0
contr_invivo <- makeContrasts(contrast = Day0-Day1, levels = colnames(coef(fit_invivo)))
fit_invivo <- contrasts.fit(fit_invivo, contr_invivo)
fit_invivo <- eBayes(fit_invivo)
tt_invivo <- topTable(fit_invivo)
tt_invivo
FDR_gsva_invivo_D1 <- tt_invivo[, c('P.Value',"adj.P.Val")]
FDR_gsva_invivo_D1$rn <- rownames(FDR_gsva_invivo_D1)
colnames(FDR_gsva_invivo_D1) <- c('P.Value.D1',"adj.P.Val.D1", 'rn')
head(FDR_gsva_invivo_D1)
# Day3 - Day0
contr_invivo <- makeContrasts(contrast = Day0-Day3, levels = colnames(coef(fit_invivo)))
fit_invivo <- contrasts.fit(fit_invivo, contr_invivo)
fit_invivo <- eBayes(fit_invivo)
tt_invivo <- topTable(fit_invivo)
tt_invivo
FDR_gsva_invivo_D3 <- tt_invivo[,c('P.Value',"adj.P.Val")]
FDR_gsva_invivo_D3$rn <- rownames(FDR_gsva_invivo_D3)
colnames(FDR_gsva_invivo_D3) <- c('P.Value.D3',"adj.P.Val.D3", 'rn')
# Day7 - Day0
contr_invivo <- makeContrasts(contrast = Day0-Day7, levels = colnames(coef(fit_invivo)))
fit_invivo <- contrasts.fit(fit_invivo, contr_invivo)
fit_invivo <- eBayes(fit_invivo)
tt_invivo <- topTable(fit_invivo)
tt_invivo
FDR_gsva_invivo_D7 <- tt_invivo[, c('P.Value',"adj.P.Val")]
FDR_gsva_invivo_D7$rn <- rownames(FDR_gsva_invivo_D7)
colnames(FDR_gsva_invivo_D7) <- c('P.Value.D7',"adj.P.Val.D7", 'rn')
# combine FDRs
FDR_gsva_invivo <- Reduce(Multimerge, list(FDR_gsva_invivo_D1, FDR_gsva_invivo_D3, FDR_gsva_invivo_D7))
head(FDR_gsva_invivo)
rownames(FDR_gsva_invivo) <- FDR_gsva_invivo$rn
FDR_gsva_invivo <- FDR_gsva_invivo[,-c(1, 2, 4, 6)]
head(FDR_gsva_invivo)


#3D
mod_3D <- model.matrix(~ factor(c(rep('Day0', 3), rep('Day4', 3), rep('Day7', 3), rep('Day10', 3), rep('Day14', 3))))
colnames(mod_3D) <- c('Day0', 'Day4', 'Day7', 'Day10', 'Day14')
fit_3D <- lmFit(gsva_3D, mod_3D)
head(coef(fit_3D))
# Day4 - Day0
contr_3D <- makeContrasts(contrast = Day0-Day4, levels = colnames(coef(fit_3D)))
fit_3D <- contrasts.fit(fit_3D, contr_3D)
fit_3D <- eBayes(fit_3D)
tt_3D <- topTable(fit_3D)
tt_3D
FDR_gsva_3D_D4 <- tt_3D[, c('P.Value',"adj.P.Val")]
FDR_gsva_3D_D4$rn <- rownames(FDR_gsva_3D_D4)
colnames(FDR_gsva_3D_D4) <- c('P.Value.D4',"adj.P.Val.D4", 'rn')
head(FDR_gsva_3D_D4)
# Day7 - Day0
contr_3D <- makeContrasts(contrast = Day0-Day7, levels = colnames(coef(fit_3D)))
fit_3D <- contrasts.fit(fit_3D, contr_3D)
fit_3D <- eBayes(fit_3D)
tt_3D <- topTable(fit_3D)
tt_3D
FDR_gsva_3D_D7 <- tt_3D[, c('P.Value',"adj.P.Val")]
FDR_gsva_3D_D7$rn <- rownames(FDR_gsva_3D_D7)
colnames(FDR_gsva_3D_D7) <- c('P.Value.D7',"adj.P.Val.D7", 'rn')
# Day10 - Day0
contr_3D <- makeContrasts(contrast = Day0-Day10, levels = colnames(coef(fit_3D)))
fit_3D <- contrasts.fit(fit_3D, contr_3D)
fit_3D <- eBayes(fit_3D)
tt_3D <- topTable(fit_3D)
tt_3D
FDR_gsva_3D_D10 <- tt_3D[, c('P.Value',"adj.P.Val")]
FDR_gsva_3D_D10$rn <- rownames(FDR_gsva_3D_D10)
colnames(FDR_gsva_3D_D10) <- c('P.Value.D10',"adj.P.Val.D10", 'rn')

# Day14 - Day0
contr_3D <- makeContrasts(contrast = Day0-Day14, levels = colnames(coef(fit_3D)))
fit_3D <- contrasts.fit(fit_3D, contr_3D)
fit_3D <- eBayes(fit_3D)
tt_3D <- topTable(fit_3D)
tt_3D
FDR_gsva_3D_D14 <- tt_3D[, c('P.Value',"adj.P.Val")]
FDR_gsva_3D_D14$rn <- rownames(FDR_gsva_3D_D14)
colnames(FDR_gsva_3D_D14) <- c('P.Value.D14',"adj.P.Val.D14", 'rn')

# combine FDRs
FDR_gsva_3D <- Reduce(Multimerge, list(FDR_gsva_3D_D4, FDR_gsva_3D_D7, FDR_gsva_3D_D10, FDR_gsva_3D_D14))
head(FDR_gsva_3D)
rownames(FDR_gsva_3D) <- FDR_gsva_3D$rn
FDR_gsva_3D <- FDR_gsva_3D[,-c(1, 2, 4, 6, 8)]
head(FDR_gsva_3D)

# all FDRs in the list
FDR_gsva_total <- list('solv in vitro' = FDR_gsva_solv, 
                       'compound in vitro' = FDR_gsva_compound, 
                       'invivo' = FDR_gsva_invivo, 
                       '3D' = FDR_gsva_3D)
FDR_gsva_total

# ggplot2 line plot with error bars and significant stars =====
library(reshape2)
library(Rmisc)
gsva_solv <- as.data.frame(t(gsva_solv))  
gsva_solv$group <- c(rep('Day0', 4), rep('Day1', 4), rep('Day2', 4), rep('Day4', 4))
gsva_solv <- melt(gsva_solv, id.vars = 'group')
head(gsva_solv)
tail(gsva_solv)
SE_gsva_solv <- summarySE(gsva_solv, measurevar = 'value', groupvars = c('group', 'variable'))
colnames(SE_gsva_solv) <- c('timepoint', 'variable', 'N', 'GSVA score', 'sd', 'se', 'ci')
SE_gsva_solv

gsva_compound <- as.data.frame(t(gsva_compound))
gsva_compound$group <- c(rep('Day0', 4), rep('Day1', 4), rep('Day2', 4), rep('Day4', 4))
gsva_compound <- melt(gsva_compound, id.vars = 'group')
head(gsva_compound)
tail(gsva_compound)
SE_gsva_compound <- summarySE(gsva_compound, measurevar = 'value', groupvars = c('group', 'variable'))
colnames(SE_gsva_compound) <- c('timepoint', 'variable', 'N', 'GSVA score', 'sd', 'se', 'ci')
SE_gsva_compound

gsva_invivo <- as.data.frame(t(gsva_invivo))
gsva_invivo$group <- c(rep('Day0', 4), rep('Day1', 2), rep('Day3', 4), rep('Day7', 4))
gsva_invivo <- melt(gsva_invivo, id.vars = 'group')
head(gsva_invivo)
tail(gsva_invivo)
SE_gsva_invivo <- summarySE(gsva_invivo, measurevar = 'value', groupvars = c('group', 'variable'))
colnames(SE_gsva_invivo) <- c('timepoint', 'variable', 'N', 'GSVA score', 'sd', 'se', 'ci')
SE_gsva_invivo

gsva_3D <- as.data.frame(t(gsva_3D))
gsva_3D$group <- c(rep('Day0', 3), rep('Day4', 3), rep('Day7', 3), rep('Day10', 3), rep('Day14', 3))
gsva_3D <- melt(gsva_3D, id.vars = 'group')
head(gsva_3D)
tail(gsva_3D)
SE_gsva_3D <- summarySE(gsva_3D, measurevar = 'value', groupvars = c('group', 'variable'))
colnames(SE_gsva_3D) <- c('timepoint', 'variable', 'N', 'GSVA score', 'sd', 'se', 'ci')
SE_gsva_3D$timepoint <- factor(SE_gsva_3D$timepoint, levels = c('Day0', 'Day4', 'Day7', 'Day10', 'Day14'))
SE_gsva_3D

# Change to signature
gsva_solv_healthy <- gsva_solv[gsva_solv$variable == 'Healthy Cell signature', ]
gsva_solv_damage <- gsva_solv[gsva_solv$variable == 'Damaged Cell signature', ]
gsva_solv_recovery <- gsva_solv[gsva_solv$variable == 'Recovery Cell signature', ]
gsva_solv_healthy$condition <- c(rep('control', 4), rep('control', 4), 
                             rep('control', 4), rep('control', 4))
gsva_solv_damage$condition <- c(rep('control', 4), rep('control', 4), 
                                rep('control', 4), rep('control', 4))
gsva_solv_recovery$condition <- c(rep('control', 4), rep('control', 4), 
                                  rep('control', 4), rep('control', 4))

gsva_compound_healthy <- gsva_compound[gsva_compound$variable == 'Healthy Cell signature', ]
gsva_compound_damage <- gsva_compound[gsva_compound$variable == 'Damaged Cell signature', ]
gsva_compound_recovery <- gsva_compound[gsva_compound$variable == 'Recovery Cell signature', ]
gsva_compound_healthy$condition <- c(rep('compound', 4), rep('compound', 4), 
                             rep('compound', 4), rep('compound', 4))
gsva_compound_damage$condition <- c(rep('compound', 4), rep('compound', 4), 
                                rep('compound', 4), rep('compound', 4))
gsva_compound_recovery$condition <- c(rep('compound', 4), rep('compound', 4), 
                                  rep('compound', 4), rep('compound', 4))

gsva_invivo_healthy <- gsva_invivo[gsva_invivo$variable == 'Healthy Cell signature', ]
gsva_invivo_damage <- gsva_invivo[gsva_invivo$variable == 'Damaged Cell signature', ]
gsva_invivo_recovery <- gsva_invivo[gsva_invivo$variable == 'Recovery Cell signature', ]
gsva_invivo_healthy$condition <- c(rep('invivo', 4), rep('invivo', 2), 
                             rep('invivo', 4), rep('invivo', 4))
gsva_invivo_damage$condition <- c(rep('invivo', 4), rep('invivo', 2), 
                        rep('invivo', 4), rep('invivo', 4))
gsva_invivo_recovery$condition <- c(rep('invivo', 4), rep('invivo', 2), 
                          rep('invivo', 4), rep('invivo', 4))


gsva_healthy <-  rbind(gsva_solv_healthy, gsva_compound_healthy, gsva_invivo_healthy)
head(gsva_healthy)
SE_gsva_healthy <- summarySE(gsva_healthy, measurevar = 'value', groupvars = c('group', 'condition'))
colnames(SE_gsva_healthy) <- c('timepoint', 'condition', 'N', 'GSVA score', 'sd', 'se', 'ci')
SE_gsva_healthy$condition <- factor(SE_gsva_healthy$condition, levels = c('invivo', 'control', 'compound'))
SE_gsva_healthy <- SE_gsva_healthy[order(SE_gsva_healthy$condition), ]
SE_gsva_healthy

gsva_damage <-  rbind(gsva_solv_damage, gsva_compound_damage, gsva_invivo_damage)
head(gsva_damage)
SE_gsva_damage <- summarySE(gsva_damage, measurevar = 'value', groupvars = c('group', 'condition'))
colnames(SE_gsva_damage) <- c('timepoint', 'condition', 'N', 'GSVA score', 'sd', 'se', 'ci')
SE_gsva_damage$condition <- factor(SE_gsva_damage$condition, levels = c('invivo', 'control', 'compound'))
SE_gsva_damage <- SE_gsva_damage[order(SE_gsva_damage$condition), ]
SE_gsva_damage

gsva_recovery <-  rbind(gsva_solv_recovery, gsva_compound_recovery, gsva_invivo_recovery)
head(gsva_recovery)
SE_gsva_recovery <- summarySE(gsva_recovery, measurevar = 'value', groupvars = c('group', 'condition'))
colnames(SE_gsva_recovery) <- c('timepoint', 'condition', 'N', 'GSVA score', 'sd', 'se', 'ci')
SE_gsva_recovery$condition <- factor(SE_gsva_recovery$condition, levels = c('invivo', 'control', 'compound'))
SE_gsva_recovery <- SE_gsva_recovery[order(SE_gsva_recovery$condition), ]
SE_gsva_recovery


library(ggplot2)
# The hex color code for the red in the plot is #F8766D.
# The hex color code for the green in the plot is #00BA38.
# The hex color code for the blue in the plot is #619CFF.
# pd <- position_dodge(0.1) # avoid errorbars overlapping
ht_gsva_healthy <- ggplot(SE_gsva_healthy, aes(timepoint, `GSVA score`, group = `condition`, color = `condition`)) + 
  geom_errorbar(aes(ymin = `GSVA score`-se, ymax = `GSVA score`+se), width = .5, position = pd) + 
  geom_line(position = pd, size=1) + geom_point(position = pd, size = 3) + 
  scale_y_continuous(limits = c(-1.0, 1.0)) + 
  # geom_text(data = data.frame(x = c('Day1', 'Day2', 'Day4'),
  #                             y = c(1.0, 0.95, 0.90)), size = 10, color = c('#F8766D', 'white', '#619CFF'),
  #           aes(x, y, label = c(rep('*', 9)))) +
  annotate('text', x = 'Day1', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', '#00BA38', '#619CFF'),
           size = 8) +
  annotate('text', x = 'Day2', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', '#00BA38', '#619CFF'),
           size = 8) +
  annotate('text', x = 'Day3', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', 'white', 'white'),
           size = 8) +
  annotate('text', x = 'Day4', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', '#00BA38', '#619CFF'),
           size = 8) +
  annotate('text', x = 'Day7', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', 'white', 'white'),
           size = 8) +
  geom_hline(yintercept = 0.0, linetype = 'dashed')+ ggtitle('Healthy Cell gene signature') + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 15/16
  ) + guides(col= guide_legend(title= "Condition"))
ht_gsva_healthy

ht_gsva_damage <- ggplot(SE_gsva_damage, aes(timepoint, `GSVA score`, group = `condition`, color = `condition`)) + 
  geom_errorbar(aes(ymin = `GSVA score`-se, ymax = `GSVA score`+se), width = .5, position = pd) + 
  geom_line(position = pd, size=1) + geom_point(position = pd, size = 3) + 
  scale_y_continuous(limits = c(-1.0, 1.0)) + 
  # geom_text(data = data.frame(x = c('Day1', 'Day2', 'Day4'),
  #                             y = c(1.0, 0.95, 0.90)), size = 10, color = c('#F8766D', 'white', '#619CFF'),
  #           aes(x, y, label = c(rep('*', 9)))) +
  annotate('text', x = 'Day1', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', '#00BA38', 'white'),
           size = 8) +
  annotate('text', x = 'Day2', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', '#00BA38', 'white'),
           size = 8) +
  annotate('text', x = 'Day3', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', 'white', 'white'),
           size = 8) +
  annotate('text', x = 'Day4', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', 'white', 'white'),
           size = 8) +
  annotate('text', x = 'Day7', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', 'white', 'white'),
           size = 8) +
  geom_hline(yintercept = 0.0, linetype = 'dashed')+ ggtitle('Damage Cell gene signature') + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 15/16
  ) + guides(col= guide_legend(title= "Condition"))
ht_gsva_damage


ht_gsva_recovery <- ggplot(SE_gsva_recovery, aes(timepoint, `GSVA score`, group = `condition`, color = `condition`)) + 
  geom_errorbar(aes(ymin = `GSVA score`-se, ymax = `GSVA score`+se), width = .5, position = pd) + 
  geom_line(position = pd, size=1) + geom_point(position = pd, size = 3) + 
  scale_y_continuous(limits = c(-1.0, 1.0)) + 
  # geom_text(data = data.frame(x = c('Day1', 'Day2', 'Day4'),
  #                             y = c(1.0, 0.95, 0.90)), size = 10, color = c('#F8766D', 'white', '#619CFF'),
  #           aes(x, y, label = c(rep('*', 9)))) +
  annotate('text', x = 'Day1', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', '#00BA38', '#619CFF'),
           size = 8) +
  annotate('text', x = 'Day2', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', '#00BA38', 'white'),
           size = 8) +
  annotate('text', x = 'Day3', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', 'white', 'white'),
           size = 8) +
  annotate('text', x = 'Day4', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('white', '#00BA38', 'white'),
           size = 8) +
  annotate('text', x = 'Day7', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', 'white', 'white'),
           size = 8) +
  geom_hline(yintercept = 0.0, linetype = 'dashed')+ ggtitle('Recovery Cell gene signature') + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") , 
    aspect.ratio = 15/16
  ) + guides(col= guide_legend(title= "Condition"))
ht_gsva_recovery

# arrange on one page
library(ggpubr)
ggarrange(ht_gsva_healthy, ht_gsva_damage, ht_gsva_recovery,
          nrow = 1, ncol = 3, common.legend = TRUE, legend = 'right') 
# (ht_gsva_solv + ht_gsva_compound)/
#   (ht_gsva_invivo +  ht_gsva_3D) + plot_layout(guides = "collect")
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/GSVAscore_Cellsignatures/')
ggplot2::ggsave(
  filename = 'GSVAscore_arranged_condition.tiff', 
  dpi = 300, 
  units = c('cm'), scale = 1
)



ht_gsva_solv <- ggplot(SE_gsva_solv, aes(timepoint, `GSVA score`, group = `variable`, color = `variable`)) + 
                        geom_errorbar(aes(ymin = `GSVA score`-se, ymax = `GSVA score`+se), width = .5, position = pd) + 
                        geom_line(position = pd, size=1) + geom_point(position = pd, size = 3) + 
                        scale_y_continuous(limits = c(-1.0, 1.0)) + 
                        # geom_text(data = data.frame(x = c('Day1', 'Day2', 'Day4'),
                        #                             y = c(1.0, 0.95, 0.90)), size = 10, color = c('#F8766D', 'white', '#619CFF'),
                        #           aes(x, y, label = c(rep('*', 9)))) +
                        annotate('text', x = 'Day1', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', '#00BA38', '#619CFF'),
                                 size = 8) +
                        annotate('text', x = 'Day2', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', '#00BA38', '#619CFF'),
                                 size = 8) +
                        annotate('text', x = 'Day4', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', 'white', '#619CFF'),
                                 size = 8) +
                        geom_hline(yintercept = 0.0, linetype = 'dashed')+ ggtitle('In vitro 2D culture without compound') + 
                        theme(
                          panel.border = element_blank(), 
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), 
                          panel.background = element_blank(), 
                          axis.line = element_line(color = "black") , 
                          aspect.ratio = 15/16
                        ) + guides(col= guide_legend(title= "Gene signatures"))
ht_gsva_solv 
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/GSVAscore_Cellsignatures/')
ggplot2::ggsave(
  filename = 'GSVAscore_invitro_2Dculture_solv.tiff', 
  dpi = 300, units = c('cm')
)

ht_gsva_compound <- ggplot(SE_gsva_compound, aes(timepoint, `GSVA score`, group = `variable`, color = `variable`)) + 
                        geom_errorbar(aes(ymin = `GSVA score`-se, ymax = `GSVA score`+se), width = .5, position = pd) + 
                        geom_line(position = pd, size=1) + geom_point(position = pd, size = 3) + 
                        scale_y_continuous(limits = c(-1.0, 1.0)) + 
                        # geom_text(data = data.frame(x = c('Day1', 'Day2', 'Day4'),
                        #                             y = c(1.0, 0.95, 0.90)), size = 10, color = c('#F8766D', 'white', '#619CFF'),
                        #           aes(x, y, label = c(rep('*', 9)))) +
                        annotate('text', x = 'Day1', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', 'white','#619CFF'),
                                 size = 8) +
                        annotate('text', x = 'Day2', y = c(0.9), label = '"*"', parse = TRUE, color= c('#F8766D'),
                                 size = 8) +
                        annotate('text', x = 'Day4', y = c(0.9), label = '"*"', parse = TRUE, color= c('#F8766D'),
                                 size = 8) +
                        geom_hline(yintercept = 0.0, linetype = 'dashed')+ ggtitle('In vitro 2D culture with compound') + 
                        theme(
                          panel.border = element_blank(), 
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), 
                          panel.background = element_blank(), 
                          axis.line = element_line(color = "black") , 
                          aspect.ratio = 15/16
                        ) + guides(col= guide_legend(title= "Gene signatures"))
ht_gsva_compound
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/GSVAscore_Cellsignatures/')
ggplot2::ggsave(
  filename = 'GSVAscore_invitro_2Dculture_compound.tiff', 
  dpi = 300, units = 'cm'
)


ht_gsva_invivo <- ggplot(SE_gsva_invivo, aes(timepoint, `GSVA score`, group = `variable`, color = `variable`)) + 
                          geom_errorbar(aes(ymin = `GSVA score`-se, ymax = `GSVA score`+se), width = .5, position = pd) + 
                          geom_line(position = pd, size=1) + geom_point(position = pd, size = 3) + 
                          scale_y_continuous(limits = c(-1.0, 1.0)) + 
                          # geom_text(data = data.frame(x = c('Day1', 'Day2', 'Day4'),
                          #                             y = c(1.0, 0.95, 0.90)), size = 10, color = c('#F8766D', 'white', '#619CFF'),
                          #           aes(x, y, label = c(rep('*', 9)))) +
                          annotate('text', x = 'Day1', y = c(0.9, 0.80), label = '"*"', parse = TRUE, color= c('#F8766D', '#00BA38'),
                                   size = 8) +
                          annotate('text', x = 'Day3', y = c(0.9, 0.80), label = '"*"', parse = TRUE, color= c('#F8766D', '#00BA38'),
                                   size = 8) +
                          annotate('text', x = 'Day7', y = c(0.9, 0.80, 0.7), label = '"*"', parse = TRUE, color= c('#F8766D', '#00BA38', '#619CFF'),
                                   size = 8) +
                          geom_hline(yintercept = 0.0, linetype = 'dashed')+ ggtitle('In vivo acute injury') + 
                          theme(
                            panel.border = element_blank(), 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), 
                            panel.background = element_blank(), 
                            axis.line = element_line(color = "black") , 
                            aspect.ratio = 15/16
                          )+ guides(col= guide_legend(title= "Gene signatures"))
ht_gsva_invivo
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/GSVAscore_Cellsignatures/')
ggplot2::ggsave(
  filename = 'GSVAscore_invivo_acute_injury.tiff', 
  dpi = 300, units = 'cm'
)

ht_gsva_3D <- ggplot(SE_gsva_3D, aes(timepoint, `GSVA score`, group = `variable`, color = `variable`)) + 
                            geom_errorbar(aes(ymin = `GSVA score`-se, ymax = `GSVA score`+se), width = .5, position = pd) + 
                            geom_line(position = pd, size=1) + geom_point(position = pd, size = 3) + 
                            scale_y_continuous(limits = c(-1.0, 1.0)) + 
                            # geom_text(data = data.frame(x = c('Day1', 'Day2', 'Day4'),
                            #                             y = c(1.0, 0.95, 0.90)), size = 10, color = c('#F8766D', 'white', '#619CFF'),
                            #           aes(x, y, label = c(rep('*', 9)))) +
                            annotate('text', x = 'Day4', y = c(0.9), label = '"*"', parse = TRUE, color= c('#F8766D'),
                                     size = 8) +
                            annotate('text', x = 'Day7', y = c(0.9), label = '"*"', parse = TRUE, color= c('#F8766D'),
                                     size = 8) +
                            annotate('text', x = 'Day10', y = c(0.9), label = '"*"', parse = TRUE, color= c('#F8766D'),
                                     size = 8) +
                            annotate('text', x = 'Day14', y = c(0.9), label = '"*"', parse = TRUE, color= c('#F8766D'),
                                     size = 8) +
                            geom_hline(yintercept = 0.0, linetype = 'dashed')+ ggtitle('Mouse 3D') + 
                            theme(
                              panel.border = element_blank(), 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), 
                              panel.background = element_blank(), 
                              axis.line = element_line(color = "black") , 
                              aspect.ratio = 15/16
                            )+ guides(col= guide_legend(title= "Gene signatures"))
ht_gsva_3D
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/GSVAscore_Cellsignatures/')
ggplot2::ggsave(
  filename = 'GSVAscore_Mouse_3D.tiff', 
  dpi = 300, units = 'cm'
)

# arrange on one page
library(ggpubr)
ggarrange(ht_gsva_solv, ht_gsva_compound, ht_gsva_invivo,ht_gsva_3D, 
          nrow = 2, ncol = 2, common.legend = TRUE, legend = 'right') 
# (ht_gsva_solv + ht_gsva_compound)/
#   (ht_gsva_invivo +  ht_gsva_3D) + plot_layout(guides = "collect")
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/GSVAscore_Cellsignatures/')
ggplot2::ggsave(
  filename = 'GSVAscore_arranged.tiff', 
  dpi = 300, 
  units = c('cm'), scale = 1
)


# Acute injury program for GSEA -------
library(DESeq2)
DEgenelist_invivo_acute_D1_up <- rownames(res_invivo_mCell_timepoint_acute[res_invivo_mCell_timepoint_acute$padj < 0.05 & 
                                                                                       res_invivo_mCell_timepoint_acute$log2FoldChange >= 1, ])
DEgenelist_invivo_acute_D1_down <- rownames(res_invivo_mCell_timepoint_acute[res_invivo_mCell_timepoint_acute$padj < 0.05 & 
                                                                             res_invivo_mCell_timepoint_acute$log2FoldChange <= -1, ])
DEgenelist_invivo_acute_D7_up <- rownames(res_invivo_mCell_timepoint_recover[res_invivo_mCell_timepoint_recover$padj < 0.05 & 
                                                                               res_invivo_mCell_timepoint_recover$log2FoldChange >= 1, ])
DEgenelist_invivo_acute_D7_down <- rownames(res_invivo_mCell_timepoint_recover[res_invivo_mCell_timepoint_recover$padj < 0.05 & 
                                                                               res_invivo_mCell_timepoint_recover$log2FoldChange <= -1, ])

length(DEgenelist_invivo_acute_D1_up)
length(DEgenelist_invivo_acute_D1_down)
length(DEgenelist_invivo_acute_D7_up)
length(DEgenelist_invivo_acute_D7_down)

# venn diagram of the overlap genes 
library(ggvenn)
DEgeneset_invivo_damaged_recovery_up <- list(`Acute injury Damaged`=DEgenelist_invivo_acute_D1_up, 
                                              `Acute injury Recovery`=DEgenelist_invivo_acute_D7_up
                                              )
ggvenn(DEgeneset_invivo_damaged_recovery_up, c('Acute injury Damaged', 'Acute injury Recovery'), auto_scale = TRUE, 
       stroke_color = 'white', stroke_size = 0.5, fill_color = c('#619CFF', '#F8766D'))
# save high resolution images
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/VennDiagram_DefinitionofAcuteInjuryProgram/')
ggplot2::ggsave(
  filename = 'venn_acuteinjury_up.tiff', 
  dpi = 300, 
  units = c('cm')
)

DEgeneset_invivo_damaged_recovery_down <- list(`Acute injury Damaged`=DEgenelist_invivo_acute_D1_down, 
                                             `Acute injury Recovery`=DEgenelist_invivo_acute_D7_down
)
ggvenn(DEgeneset_invivo_damaged_recovery_down, c('Acute injury Damaged', 'Acute injury Recovery'), auto_scale = TRUE, 
       stroke_color = 'white', stroke_size = 0.5, fill_color = c('#619CFF', '#F8766D')) 

# save high resolution imagesm 
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/VennDiagram_DefinitionofAcuteInjuryProgram/')
ggplot2::ggsave(
  filename = 'venn_acuteinjury_down.tiff', 
  dpi = 300, 
  units = c('cm')
)


# Extract the overlap genes 
overlap_DEgene_invivo_damaged_recovery_up <- Reduce(intersect, list(DEgenelist_invivo_acute_D1_up, DEgenelist_invivo_acute_D7_up))
overlap_DEgene_invivo_damaged_recovery_down <- Reduce(intersect, list(DEgenelist_invivo_acute_D1_down, DEgenelist_invivo_acute_D7_down))

# heatmap of overlap genes
library(ComplexHeatmap)
library(circlize)
overlap_DEgene_invivo_damaged_recovery_up <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                        overlap_DEgene_invivo_damaged_recovery_up, ]
nrow(overlap_DEgene_invivo_damaged_recovery_up)
head(overlap_DEgene_invivo_damaged_recovery_up)
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/GSEA_chronicinjury/')
write.csv(overlap_DEgene_invivo_damaged_recovery_up, file = 'overlap_DEgene_invivo_damaged_recovery_up.csv')

overlap_DEgene_invivo_damaged_recovery_down <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                 overlap_DEgene_invivo_damaged_recovery_down, ]
nrow(overlap_DEgene_invivo_damaged_recovery_down)
write.csv(overlap_DEgene_invivo_damaged_recovery_down, file = 'overlap_DEgene_invivo_damaged_recovery_down.csv')

overlap_DEgene_invivo_damaged_recovery_up <- as.data.frame(t(scale(t(overlap_DEgene_invivo_damaged_recovery_up))))
head(overlap_DEgene_invivo_damaged_recovery_up)

overlap_DEgene_invivo_damaged_recovery_up_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_up, 
                                                         `D0` = apply(overlap_DEgene_invivo_damaged_recovery_up[1:4], 1, mean))
overlap_DEgene_invivo_damaged_recovery_up_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_up_mean, 
                                                         `D1` = apply(overlap_DEgene_invivo_damaged_recovery_up[5: 6], 1, mean))
overlap_DEgene_invivo_damaged_recovery_up_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_up_mean, 
                                                         `D3` = apply(overlap_DEgene_invivo_damaged_recovery_up[7: 10], 1, mean))
overlap_DEgene_invivo_damaged_recovery_up_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_up_mean, 
                                                         `D7` = apply(overlap_DEgene_invivo_damaged_recovery_up[11: 14], 1, mean))

overlap_DEgene_invivo_damaged_recovery_up_mean <- overlap_DEgene_invivo_damaged_recovery_up_mean[,c('D0', 'D1', 'D3', 'D7')]
head(overlap_DEgene_invivo_damaged_recovery_up_mean)

ht_overlap_DEgene_invivo_damaged_recovery_up <- Heatmap(overlap_DEgene_invivo_damaged_recovery_up_mean,
                                                        col = color_ht,
                                                        width = unit(8*ncol(overlap_DEgene_invivo_damaged_recovery_up_mean), 'mm'),
                                                        row_title = 'Recovery gene signature', 
                                                        cluster_columns = FALSE,
                                                        cluster_rows = FALSE, 
                                                        # cluster_rows = FALSE,
                                                        show_row_names = FALSE,
                                                        column_names_side = 'top', 
                                                        # rect_gp = gpar(col = 'white', lwd = 1), 
                                                        show_heatmap_legend = FALSE
                                                        
                              )
ht_overlap_DEgene_invivo_damaged_recovery_up

overlap_DEgene_invivo_damaged_recovery_down <- as.data.frame(t(scale(t(overlap_DEgene_invivo_damaged_recovery_down))))
head(overlap_DEgene_invivo_damaged_recovery_down)

overlap_DEgene_invivo_damaged_recovery_down_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_down, 
                                                         `D0` = apply(overlap_DEgene_invivo_damaged_recovery_down[1:4], 1, mean))
overlap_DEgene_invivo_damaged_recovery_down_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_down_mean, 
                                                         `D1` = apply(overlap_DEgene_invivo_damaged_recovery_down[5: 6], 1, mean))
overlap_DEgene_invivo_damaged_recovery_down_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_down_mean, 
                                                         `D3` = apply(overlap_DEgene_invivo_damaged_recovery_down[7: 10], 1, mean))
overlap_DEgene_invivo_damaged_recovery_down_mean <- mutate(overlap_DEgene_invivo_damaged_recovery_down_mean, 
                                                         `D7` = apply(overlap_DEgene_invivo_damaged_recovery_down[11: 14], 1, mean))

overlap_DEgene_invivo_damaged_recovery_down_mean <- overlap_DEgene_invivo_damaged_recovery_down_mean[,c('D0', 'D1', 'D3', 'D7')]
head(overlap_DEgene_invivo_damaged_recovery_down_mean)

ht_overlap_DEgene_invivo_damaged_recovery_down <- Heatmap(overlap_DEgene_invivo_damaged_recovery_down_mean,
                                                        col = color_ht,
                                                        width = unit(8*ncol(overlap_DEgene_invivo_damaged_recovery_down_mean), 'mm'),
                                                        row_title = 'Recovery gene signature', 
                                                        cluster_columns = FALSE,
                                                        cluster_rows = FALSE, 
                                                        # cluster_rows = FALSE,
                                                        show_row_names = FALSE,
                                                        column_names_side = 'top', 
                                                        # rect_gp = gpar(col = 'white', lwd = 1), 
                                                        show_heatmap_legend = FALSE
                                                        
)
ht_overlap_DEgene_invivo_damaged_recovery_down

# GO analysis of the overlapped genes
# library(clusterProfiler)
# library(org.Mm.eg.db)
# 
# GO_overlap_DEgene_invivo_damaged_recovery_up <- enrichGO(gene = rownames(overlap_DEgene_invivo_damaged_recovery_up), 
#                                                          OrgDb = org.Mm.eg.db, 
#                                                          keyType = 'SYMBOL', 
#                                                          ont = 'BP', 
#                                                          pvalueCutoff = 0.05, 
#                                                          pAdjustMethod = 'BH',
#                                                          qvalueCutoff = 0.05,
#                                                          readable = TRUE
#                                                         )
# head(GO_overlap_DEgene_invivo_damaged_recovery_up)
# 
# GO_overlap_DEgene_invivo_damaged_recovery_down <- enrichGO(gene = rownames(overlap_DEgene_invivo_damaged_recovery_down), 
#                                                            OrgDb = org.Mm.eg.db, 
#                                                            keyType = 'SYMBOL', 
#                                                            ont = 'BP', 
#                                                            pvalueCutoff = 0.05, 
#                                                            pAdjustMethod = 'BH',
#                                                            qvalueCutoff = 0.05,
#                                                            readable = TRUE
#                                                           )
# head(GO_overlap_DEgene_invivo_damaged_recovery_down)

# dotplot(GO_overlap_DEgene_invivo_damaged_recovery_down, showCategory = 15)
# rm(GO_overlap_DEgene_invivo_damaged_recovery_down)
# rm(GO_overlap_DEgene_invivo_damaged_recovery_up)

library(xlsx)
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/GOanalysis_DAVID/')
GO_overlap_DEgene_invivo_damaged_recovery_total <- read.xlsx('overlap_invivo_damaged_recovery_total_GOBP_forR.xlsx', 1)
head(GO_overlap_DEgene_invivo_damaged_recovery_total)
GO_overlap_DEgene_invivo_damaged_recovery_total$Count <- as.numeric(GO_overlap_DEgene_invivo_damaged_recovery_total$Count)
GO_overlap_DEgene_invivo_damaged_recovery_total$GeneRatio <- as.numeric(GO_overlap_DEgene_invivo_damaged_recovery_total$GeneRatio)
GO_overlap_DEgene_invivo_damaged_recovery_total$FDR <- as.numeric(GO_overlap_DEgene_invivo_damaged_recovery_total$FDR)

library(ggplot2)
ggplot(data = GO_overlap_DEgene_invivo_damaged_recovery_total, aes(x = GeneRatio, y = factor(Description, levels = rev(Description)), 
                                                                   color = FDR, size = Count)) + 
  geom_point() + 
  scale_color_gradient(high = '#619CFF', low = '#F8766D') + 
  theme_bw() + 
  ylab('') + xlab('GeneRatio') + theme(axis.text.x = element_text(size = 12), 
                                      axis.text.y = element_text(size = 12, face = 'bold'), 
                                      axis.title.x = element_text(size = 12), 
                                      legend.title = element_text(size = 12), 
                                      legend.text = element_text(size = 12)) + 
  scale_size(range = c(5,10)) + 
  scale_x_continuous(limits = c(0.03,0.11))

# save high resolution imagesm 
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/PlotsFromR/')
ggplot2::ggsave(
  filename = 'GO_DAVID_dotplot.tiff', 
  dpi = 300, 
  units = c('cm')
)

# GSEA report plot from GSEA software -------
# library(Rtoolbox)
# replotGSEA(path = 'D:/foldername/1. mCell_Dedifferentiation/RawData/GSEA_chronicinjury/GSEA_Results/compound_D2D0.Gsea.1700490958963/', 
#            gene.set = 'ACUTE INJURY DE GENE SET', class.name = c('D2 D0'))


# chronic injury (2w, 4w from scRNA data) create GSEA pathways ------
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/')
library(xlsx)
write.xlsx(normcounts_2D_mCell_timepoint_solv_D1, file = 'normcounts_2D_mCell_solv.xlsx', col.names = TRUE, row.names = TRUE)
write.xlsx(normcounts_2D_mCell_timepoint_compound_D1, file = 'normcounts_2D_mCell_compound.xlsx', col.names = TRUE, row.names = TRUE)
write.csv(normcounts_invivo_mCell_timepoint_acute, file = 'normcounts_invivoacute.csv', row.names = TRUE)
write.csv(normcounts_3D, file = 'normcounts_3D.csv', row.names = TRUE)

# Pathways heatmap of genes in Cell capillarization with compound -----
compound_eNOS_pathway <- c('Nos3', 'Gucy1a1', 'Gucy1b1') # PMID: 25131509
Hedgehog_pathway <- c('Ptch1', 'Gli2', 'Gli3', 'Hhip', 'Twist2', 'Sfrp1', 'Spp1') #  PMID: 22362915
endocytosis_receptor <- c('Cd206', 'Stab1', 'Stab2', 'Fcgr2b') # PMID: 25131509

# compound_eNOS_pathway ====
Gene_compound_eNOS_pathway_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                           compound_eNOS_pathway, ]

head(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_invitro)

Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv <- Gene_compound_eNOS_pathway_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv)

Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound <- Gene_compound_eNOS_pathway_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound)


Gene_compound_eNOS_pathway_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                    compound_eNOS_pathway, ]
head(Gene_compound_eNOS_pathway_rld_counts_invivo)


Gene_compound_eNOS_pathway_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                        compound_eNOS_pathway, ]
# Gene_compound_eNOS_pathway_rld_counts_3D <- as.data.frame(t(Gene_compound_eNOS_pathway_rld_counts_3D))
# Gene_compound_eNOS_pathway_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_compound_eNOS_pathway_rld_counts_3D)


Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv))))
head(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv)

Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv, 
                                                                   `D0` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean, 
                                                                   `D1` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean, 
                                                                   `D2` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean, 
                                                                   `D4` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean <- Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean)

Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound))))
head(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound)

Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound, 
                                                                   `D0` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean, 
                                                                   `D1` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean, 
                                                                   `D2` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean, 
                                                                   `D4` = apply(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean <- Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean)



Gene_compound_eNOS_pathway_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_compound_eNOS_pathway_rld_counts_invivo))))
head(Gene_compound_eNOS_pathway_rld_counts_invivo)

Gene_compound_eNOS_pathway_rld_counts_invivo_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_invivo, 
                                                            `D0` = apply(Gene_compound_eNOS_pathway_rld_counts_invivo[1:4], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_invivo_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_invivo_mean, 
                                                            `D1` = apply(Gene_compound_eNOS_pathway_rld_counts_invivo[5: 6], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_invivo_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_invivo_mean, 
                                                            `D3` = apply(Gene_compound_eNOS_pathway_rld_counts_invivo[7: 10], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_invivo_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_invivo_mean, 
                                                            `D7` = apply(Gene_compound_eNOS_pathway_rld_counts_invivo[11: 14], 1, mean))

Gene_compound_eNOS_pathway_rld_counts_invivo_mean <- Gene_compound_eNOS_pathway_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_compound_eNOS_pathway_rld_counts_invivo_mean)


Gene_compound_eNOS_pathway_rld_counts_3D <- as.data.frame(t(scale(t(Gene_compound_eNOS_pathway_rld_counts_3D))))

Gene_compound_eNOS_pathway_rld_counts_3D_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_3D, 
                                                              `D0` = apply(Gene_compound_eNOS_pathway_rld_counts_3D[1: 3], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_3D_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_3D_mean, 
                                                              `D4` = apply(Gene_compound_eNOS_pathway_rld_counts_3D_mean[4: 6], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_3D_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_3D_mean, 
                                                              `D7` = apply(Gene_compound_eNOS_pathway_rld_counts_3D_mean[7: 9], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_3D_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_3D_mean, 
                                                              `D10` = apply(Gene_compound_eNOS_pathway_rld_counts_3D_mean[10: 12], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_3D_mean <- mutate(Gene_compound_eNOS_pathway_rld_counts_3D_mean, 
                                                              `D14` = apply(Gene_compound_eNOS_pathway_rld_counts_3D_mean[13: 15], 1, mean))
Gene_compound_eNOS_pathway_rld_counts_3D_mean <- Gene_compound_eNOS_pathway_rld_counts_3D_mean[,16: 20]
head(Gene_compound_eNOS_pathway_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_compound_eNOS_pathway_solv <- Heatmap(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean,
                                    col = color_ht,
                                    width = unit(8*ncol(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean), 'mm'), 
                                    height = unit(6* nrow(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean), 'mm'),
                                    column_title = 'In vitro 2D without compound',
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE, 
                                    # cluster_rows = FALSE,
                                    # show_row_names = FALSE,
                                    column_names_side = 'top', 
                                    rect_gp = gpar(col = 'white', lwd = 1), 
                                    heatmap_legend_param = list(title = 'Normalized counts')
                                    # column_split = column_split_overlap,
                                    # row_names_max_width = max_text_width(rownames(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_solv_mean,
                                    #                                      gp = gpar(fontsize = 12))    )
)
ht_compound_eNOS_pathway_solv

ht_compound_eNOS_pathway_compound <- Heatmap(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean,
                                    col = color_ht,
                                    width = unit(8*ncol(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean), 'mm'), 
                                    height = unit(6* nrow(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean), 'mm'),
                                    column_title = 'In vitro 2D with compound',
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE, 
                                    # cluster_rows = FALSE,
                                    # show_row_names = FALSE,
                                    column_names_side = 'top', 
                                    rect_gp = gpar(col = 'white', lwd = 1), 
                                    heatmap_legend_param = list(title = 'Normalized counts')
                                    # column_split = column_split_overlap,
                                    # row_names_max_width = max_text_width(rownames(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_compound_mean,
                                    #                                      gp = gpar(fontsize = 12))    )
)
ht_compound_eNOS_pathway_compound


ht_compound_eNOS_pathway_invivo <- Heatmap(Gene_compound_eNOS_pathway_rld_counts_invivo_mean,
                                      col = color_ht,
                                      width = unit(8*ncol(Gene_compound_eNOS_pathway_rld_counts_invivo_mean), 'mm'), 
                                      height = unit(6* nrow(Gene_compound_eNOS_pathway_rld_counts_invivo_mean), 'mm'),
                                      column_title = 'In vivo acute injury',
                                      cluster_columns = FALSE,
                                      cluster_rows = FALSE, 
                                      # cluster_rows = FALSE,
                                      # show_row_names = FALSE,
                                      column_names_side = 'top', 
                                      rect_gp = gpar(col = 'white', lwd = 1), 
                                      heatmap_legend_param = list(title = 'Normalized counts')
                                      # column_split = column_split_overlap,
                                      # row_names_max_width = max_text_width(rownames(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_invivo_mean,
                                      #                                      gp = gpar(fontsize = 12))    )
)
ht_compound_eNOS_pathway_invivo

ht_compound_eNOS_pathway_3D <- Heatmap(Gene_compound_eNOS_pathway_rld_counts_3D_mean,
                                        col = color_ht,
                                        width = unit(8*ncol(Gene_compound_eNOS_pathway_rld_counts_3D_mean), 'mm'), 
                                        height = unit(6* nrow(Gene_compound_eNOS_pathway_rld_counts_3D_mean), 'mm'),
                                        column_title = 'Moues 3D',
                                        cluster_columns = FALSE,
                                        cluster_rows = FALSE, 
                                        # cluster_rows = FALSE,
                                        # show_row_names = FALSE,
                                        column_names_side = 'top', 
                                        rect_gp = gpar(col = 'white', lwd = 1), 
                                        heatmap_legend_param = list(title = 'Normalized counts')
                                        # column_split = column_split_overlap,
                                        # row_names_max_width = max_text_width(rownames(Gene_compound_eNOS_pathway_rld_counts_2D_mCell_3D_mean,
                                        #                                      gp = gpar(fontsize = 12))    )
)

# Hedgehog_pathway ====
Gene_Hedgehog_pathway_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                           Hedgehog_pathway, ]

head(Gene_Hedgehog_pathway_rld_counts_2D_mCell_invitro)

Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv <- Gene_Hedgehog_pathway_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv)

Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound <- Gene_Hedgehog_pathway_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound)


Gene_Hedgehog_pathway_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                                Hedgehog_pathway, ]
head(Gene_Hedgehog_pathway_rld_counts_invivo)


Gene_Hedgehog_pathway_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                    Hedgehog_pathway, ]
# Gene_Hedgehog_pathway_rld_counts_3D <- as.data.frame(t(Gene_Hedgehog_pathway_rld_counts_3D))
# Gene_Hedgehog_pathway_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_Hedgehog_pathway_rld_counts_3D)


Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv))))
head(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv)

Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv, 
                                                               `D0` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean, 
                                                               `D1` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean, 
                                                               `D2` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean, 
                                                               `D4` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean <- Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean)

Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound))))
head(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound)

Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound, 
                                                               `D0` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean, 
                                                               `D1` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean, 
                                                               `D2` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean, 
                                                               `D4` = apply(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean <- Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean)



Gene_Hedgehog_pathway_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_Hedgehog_pathway_rld_counts_invivo))))
head(Gene_Hedgehog_pathway_rld_counts_invivo)

Gene_Hedgehog_pathway_rld_counts_invivo_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_invivo, 
                                                        `D0` = apply(Gene_Hedgehog_pathway_rld_counts_invivo[1:4], 1, mean))
Gene_Hedgehog_pathway_rld_counts_invivo_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_invivo_mean, 
                                                        `D1` = apply(Gene_Hedgehog_pathway_rld_counts_invivo[5: 6], 1, mean))
Gene_Hedgehog_pathway_rld_counts_invivo_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_invivo_mean, 
                                                        `D3` = apply(Gene_Hedgehog_pathway_rld_counts_invivo[7: 10], 1, mean))
Gene_Hedgehog_pathway_rld_counts_invivo_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_invivo_mean, 
                                                        `D7` = apply(Gene_Hedgehog_pathway_rld_counts_invivo[11: 14], 1, mean))

Gene_Hedgehog_pathway_rld_counts_invivo_mean <- Gene_Hedgehog_pathway_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_Hedgehog_pathway_rld_counts_invivo_mean)


Gene_Hedgehog_pathway_rld_counts_3D <- as.data.frame(t(scale(t(Gene_Hedgehog_pathway_rld_counts_3D))))

Gene_Hedgehog_pathway_rld_counts_3D_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_3D, 
                                                          `D0` = apply(Gene_Hedgehog_pathway_rld_counts_3D[1: 3], 1, mean))
Gene_Hedgehog_pathway_rld_counts_3D_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_3D_mean, 
                                                          `D4` = apply(Gene_Hedgehog_pathway_rld_counts_3D_mean[4: 6], 1, mean))
Gene_Hedgehog_pathway_rld_counts_3D_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_3D_mean, 
                                                          `D7` = apply(Gene_Hedgehog_pathway_rld_counts_3D_mean[7: 9], 1, mean))
Gene_Hedgehog_pathway_rld_counts_3D_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_3D_mean, 
                                                          `D10` = apply(Gene_Hedgehog_pathway_rld_counts_3D_mean[10: 12], 1, mean))
Gene_Hedgehog_pathway_rld_counts_3D_mean <- mutate(Gene_Hedgehog_pathway_rld_counts_3D_mean, 
                                                          `D14` = apply(Gene_Hedgehog_pathway_rld_counts_3D_mean[13: 15], 1, mean))
Gene_Hedgehog_pathway_rld_counts_3D_mean <- Gene_Hedgehog_pathway_rld_counts_3D_mean[,16: 20]
head(Gene_Hedgehog_pathway_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_Hedgehog_pathway_solv <- Heatmap(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean,
                                     col = color_ht,
                                     width = unit(8*ncol(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean), 'mm'), 
                                     height = unit(6* nrow(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean), 'mm'),
                                     column_title = 'In vitro 2D without compound',
                                     cluster_columns = FALSE,
                                     cluster_rows = FALSE, 
                                     # cluster_rows = FALSE,
                                     # show_row_names = FALSE,
                                     column_names_side = 'top', 
                                     rect_gp = gpar(col = 'white', lwd = 1), 
                                     heatmap_legend_param = list(title = 'Normalized counts')
                                     # column_split = column_split_overlap,
                                     # row_names_max_width = max_text_width(rownames(Gene_Hedgehog_pathway_rld_counts_2D_mCell_solv_mean,
                                     #                                      gp = gpar(fontsize = 12))    )
)
ht_Hedgehog_pathway_solv

ht_Hedgehog_pathway_compound <- Heatmap(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean,
                                     col = color_ht,
                                     width = unit(8*ncol(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean), 'mm'), 
                                     height = unit(6* nrow(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean), 'mm'),
                                     column_title = 'In vitro 2D with compound',
                                     cluster_columns = FALSE,
                                     cluster_rows = FALSE, 
                                     # cluster_rows = FALSE,
                                     # show_row_names = FALSE,
                                     column_names_side = 'top', 
                                     rect_gp = gpar(col = 'white', lwd = 1), 
                                     heatmap_legend_param = list(title = 'Normalized counts')
                                     # column_split = column_split_overlap,
                                     # row_names_max_width = max_text_width(rownames(Gene_Hedgehog_pathway_rld_counts_2D_mCell_compound_mean,
                                     #                                      gp = gpar(fontsize = 12))    )
)
ht_Hedgehog_pathway_compound


ht_Hedgehog_pathway_invivo <- Heatmap(Gene_Hedgehog_pathway_rld_counts_invivo_mean,
                                       col = color_ht,
                                       width = unit(8*ncol(Gene_Hedgehog_pathway_rld_counts_invivo_mean), 'mm'), 
                                       height = unit(6* nrow(Gene_Hedgehog_pathway_rld_counts_invivo_mean), 'mm'),
                                       column_title = 'In vivo acute injury',
                                       cluster_columns = FALSE,
                                       cluster_rows = FALSE, 
                                       # cluster_rows = FALSE,
                                       # show_row_names = FALSE,
                                       column_names_side = 'top', 
                                       rect_gp = gpar(col = 'white', lwd = 1), 
                                       heatmap_legend_param = list(title = 'Normalized counts')
                                       # column_split = column_split_overlap,
                                       # row_names_max_width = max_text_width(rownames(Gene_Hedgehog_pathway_rld_counts_2D_mCell_invivo_mean,
                                       #                                      gp = gpar(fontsize = 12))    )
)
ht_Hedgehog_pathway_invivo

ht_Hedgehog_pathway_3D <- Heatmap(Gene_Hedgehog_pathway_rld_counts_3D_mean,
                                         col = color_ht,
                                         width = unit(8*ncol(Gene_Hedgehog_pathway_rld_counts_3D_mean), 'mm'), 
                                         height = unit(6* nrow(Gene_Hedgehog_pathway_rld_counts_3D_mean), 'mm'),
                                         column_title = 'Moues 3D',
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE, 
                                         # cluster_rows = FALSE,
                                         # show_row_names = FALSE,
                                         column_names_side = 'top', 
                                         rect_gp = gpar(col = 'white', lwd = 1), 
                                         heatmap_legend_param = list(title = 'Normalized counts')
                                         # column_split = column_split_overlap,
                                         # row_names_max_width = max_text_width(rownames(Gene_Hedgehog_pathway_rld_counts_2D_mCell_3D_mean,
                                         #                                      gp = gpar(fontsize = 12))    )
)

# endocytosis_receptor ====
Gene_endocytosis_receptor_rld_counts_2D_mCell_invitro <- rld_counts_invitro[rownames(rld_counts_invitro) %in% 
                                                                          endocytosis_receptor, ]

head(Gene_endocytosis_receptor_rld_counts_2D_mCell_invitro)

Gene_endocytosis_receptor_rld_counts_2D_mCell_solv <- Gene_endocytosis_receptor_rld_counts_2D_mCell_invitro[,1:16]
head(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv)

Gene_endocytosis_receptor_rld_counts_2D_mCell_compound <- Gene_endocytosis_receptor_rld_counts_2D_mCell_invitro[,c(1:4, 17:28)]
head(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound)


Gene_endocytosis_receptor_rld_counts_invivo <- rld_counts_invivo[rownames(rld_counts_invivo) %in% 
                                                               endocytosis_receptor, ]
head(Gene_endocytosis_receptor_rld_counts_invivo)


Gene_endocytosis_receptor_rld_counts_3D <- rld_counts_3D[rownames(rld_counts_3D) %in% 
                                                                   endocytosis_receptor, ]
# Gene_endocytosis_receptor_rld_counts_3D <- as.data.frame(t(Gene_endocytosis_receptor_rld_counts_3D))
# Gene_endocytosis_receptor_rld_counts_3D$group <- c(rep('0', 3), rep('4', 3), rep('7', 3), rep('10', 3), rep('14', 3))
head(Gene_endocytosis_receptor_rld_counts_3D)


Gene_endocytosis_receptor_rld_counts_2D_mCell_solv <- as.data.frame(t(scale(t(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv))))
head(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv)

Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv, 
                                                              `D0` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv[1:4], 1, mean))
Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean, 
                                                              `D1` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv[5: 8], 1, mean))
Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean, 
                                                              `D2` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv[9: 12], 1, mean))
Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean, 
                                                              `D4` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv[13: 16], 1, mean))

Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean <- Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean)

Gene_endocytosis_receptor_rld_counts_2D_mCell_compound <- as.data.frame(t(scale(t(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound))))
head(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound)

Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound, 
                                                              `D0` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound[1:4], 1, mean))
Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean, 
                                                              `D1` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound[5: 8], 1, mean))
Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean, 
                                                              `D2` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound[9: 12], 1, mean))
Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean <- mutate(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean, 
                                                              `D4` = apply(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound[13: 16], 1, mean))

Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean <- Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean[,c('D0', 'D1', 'D2', 'D4')]
head(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean)



Gene_endocytosis_receptor_rld_counts_invivo <- as.data.frame(t(scale(t(Gene_endocytosis_receptor_rld_counts_invivo))))
head(Gene_endocytosis_receptor_rld_counts_invivo)

Gene_endocytosis_receptor_rld_counts_invivo_mean <- mutate(Gene_endocytosis_receptor_rld_counts_invivo, 
                                                       `D0` = apply(Gene_endocytosis_receptor_rld_counts_invivo[1:4], 1, mean))
Gene_endocytosis_receptor_rld_counts_invivo_mean <- mutate(Gene_endocytosis_receptor_rld_counts_invivo_mean, 
                                                       `D1` = apply(Gene_endocytosis_receptor_rld_counts_invivo[5: 6], 1, mean))
Gene_endocytosis_receptor_rld_counts_invivo_mean <- mutate(Gene_endocytosis_receptor_rld_counts_invivo_mean, 
                                                       `D3` = apply(Gene_endocytosis_receptor_rld_counts_invivo[7: 10], 1, mean))
Gene_endocytosis_receptor_rld_counts_invivo_mean <- mutate(Gene_endocytosis_receptor_rld_counts_invivo_mean, 
                                                       `D7` = apply(Gene_endocytosis_receptor_rld_counts_invivo[11: 14], 1, mean))

Gene_endocytosis_receptor_rld_counts_invivo_mean <- Gene_endocytosis_receptor_rld_counts_invivo_mean[,c('D0', 'D1', 'D3', 'D7')]
head(Gene_endocytosis_receptor_rld_counts_invivo_mean)


Gene_endocytosis_receptor_rld_counts_3D <- as.data.frame(t(scale(t(Gene_endocytosis_receptor_rld_counts_3D))))

Gene_endocytosis_receptor_rld_counts_3D_mean <- mutate(Gene_endocytosis_receptor_rld_counts_3D, 
                                                         `D0` = apply(Gene_endocytosis_receptor_rld_counts_3D[1: 3], 1, mean))
Gene_endocytosis_receptor_rld_counts_3D_mean <- mutate(Gene_endocytosis_receptor_rld_counts_3D_mean, 
                                                         `D4` = apply(Gene_endocytosis_receptor_rld_counts_3D_mean[4: 6], 1, mean))
Gene_endocytosis_receptor_rld_counts_3D_mean <- mutate(Gene_endocytosis_receptor_rld_counts_3D_mean, 
                                                         `D7` = apply(Gene_endocytosis_receptor_rld_counts_3D_mean[7: 9], 1, mean))
Gene_endocytosis_receptor_rld_counts_3D_mean <- mutate(Gene_endocytosis_receptor_rld_counts_3D_mean, 
                                                         `D10` = apply(Gene_endocytosis_receptor_rld_counts_3D_mean[10: 12], 1, mean))
Gene_endocytosis_receptor_rld_counts_3D_mean <- mutate(Gene_endocytosis_receptor_rld_counts_3D_mean, 
                                                         `D14` = apply(Gene_endocytosis_receptor_rld_counts_3D_mean[13: 15], 1, mean))
Gene_endocytosis_receptor_rld_counts_3D_mean <- Gene_endocytosis_receptor_rld_counts_3D_mean[,16: 20]
head(Gene_endocytosis_receptor_rld_counts_3D_mean)

library(ComplexHeatmap)
library(circlize)
ht_endocytosis_receptor_solv <- Heatmap(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean,
                                    col = color_ht,
                                    width = unit(8*ncol(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean), 'mm'), 
                                    height = unit(6* nrow(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean), 'mm'),
                                    column_title = 'In vitro 2D without compound',
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE, 
                                    # cluster_rows = FALSE,
                                    # show_row_names = FALSE,
                                    column_names_side = 'top', 
                                    rect_gp = gpar(col = 'white', lwd = 1), 
                                    heatmap_legend_param = list(title = 'Normalized counts')
                                    # column_split = column_split_overlap,
                                    # row_names_max_width = max_text_width(rownames(Gene_endocytosis_receptor_rld_counts_2D_mCell_solv_mean,
                                    #                                      gp = gpar(fontsize = 12))    )
)
ht_endocytosis_receptor_solv

ht_endocytosis_receptor_compound <- Heatmap(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean,
                                    col = color_ht,
                                    width = unit(8*ncol(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean), 'mm'), 
                                    height = unit(6* nrow(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean), 'mm'),
                                    column_title = 'In vitro 2D with compound',
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE, 
                                    # cluster_rows = FALSE,
                                    # show_row_names = FALSE,
                                    column_names_side = 'top', 
                                    rect_gp = gpar(col = 'white', lwd = 1), 
                                    heatmap_legend_param = list(title = 'Normalized counts')
                                    # column_split = column_split_overlap,
                                    # row_names_max_width = max_text_width(rownames(Gene_endocytosis_receptor_rld_counts_2D_mCell_compound_mean,
                                    #                                      gp = gpar(fontsize = 12))    )
)
ht_endocytosis_receptor_compound


ht_endocytosis_receptor_invivo <- Heatmap(Gene_endocytosis_receptor_rld_counts_invivo_mean,
                                      col = color_ht,
                                      width = unit(8*ncol(Gene_endocytosis_receptor_rld_counts_invivo_mean), 'mm'), 
                                      height = unit(6* nrow(Gene_endocytosis_receptor_rld_counts_invivo_mean), 'mm'),
                                      column_title = 'In vivo acute injury',
                                      cluster_columns = FALSE,
                                      cluster_rows = FALSE, 
                                      # cluster_rows = FALSE,
                                      # show_row_names = FALSE,
                                      column_names_side = 'top', 
                                      rect_gp = gpar(col = 'white', lwd = 1), 
                                      heatmap_legend_param = list(title = 'Normalized counts')
                                      # column_split = column_split_overlap,
                                      # row_names_max_width = max_text_width(rownames(Gene_endocytosis_receptor_rld_counts_2D_mCell_invivo_mean,
                                      #                                      gp = gpar(fontsize = 12))    )
)
ht_endocytosis_receptor_invivo

ht_endocytosis_receptor_3D <- Heatmap(Gene_endocytosis_receptor_rld_counts_3D_mean,
                                        col = color_ht,
                                        width = unit(8*ncol(Gene_endocytosis_receptor_rld_counts_3D_mean), 'mm'), 
                                        height = unit(6* nrow(Gene_endocytosis_receptor_rld_counts_3D_mean), 'mm'),
                                        column_title = 'Moues 3D',
                                        cluster_columns = FALSE,
                                        cluster_rows = FALSE, 
                                        # cluster_rows = FALSE,
                                        # show_row_names = FALSE,
                                        column_names_side = 'top', 
                                        rect_gp = gpar(col = 'white', lwd = 1), 
                                        heatmap_legend_param = list(title = 'Normalized counts')
                                        # column_split = column_split_overlap,
                                        # row_names_max_width = max_text_width(rownames(Gene_endocytosis_receptor_rld_counts_2D_mCell_3D_mean,
                                        #                                      gp = gpar(fontsize = 12))    )
)

# combine ht =====
draw(ht_compound_eNOS_pathway_solv %v% ht_Hedgehog_pathway_solv %v% ht_endocytosis_receptor_solv)
draw(ht_compound_eNOS_pathway_compound %v% ht_Hedgehog_pathway_compound %v% ht_endocytosis_receptor_compound)
draw(ht_compound_eNOS_pathway_invivo %v% ht_Hedgehog_pathway_invivo %v% ht_endocytosis_receptor_invivo)
draw(ht_compound_eNOS_pathway_3D %v% ht_Hedgehog_pathway_3D %v% ht_endocytosis_receptor_3D)

# Pathway heatmap - Tsuchida 2018 -----
library(ComplexHeatmap)
library(circlize)

# GSEA in vivo acute injury D1 vs Healthy
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/GSEA/GSEA_Results/20231208_invivo_D1.Gsea.1702031036775/')

GSEA_invivo_D1_up <- read.table('gsea_report_for_D1_1702031036775.tsv', header = TRUE, sep = '\t')
GSEA_invivo_D1_up$Regulation <- rep(c('UP'), 31)
GSEA_invivo_D1_down <- read.table('gsea_report_for_D0_1702031036775.tsv', header = TRUE, sep = '\t')
nrow(GSEA_invivo_D1_down)
GSEA_invivo_D1_down$Regulation <- rep(c('DOWN'), 19)
GSEA_invivo_D1 <- rbind(GSEA_invivo_D1_up, GSEA_invivo_D1_down)
GSEA_invivo_D1 <- GSEA_invivo_D1[, colnames(GSEA_invivo_D1) %in% c('NAME', 'SIZE', 'ES', 'NES', 'NOM.p.val', 'FDR.q.val', 
                                                                   'Regulation')]
head(GSEA_invivo_D1)
tail(GSEA_invivo_D1)
GSEA_invivo_D1_sig <- GSEA_invivo_D1[GSEA_invivo_D1$FDR.q.val < 0.01, ]
GSEA_invivo_D1_sig$NAME <- gsub('HALLMARK_', '', GSEA_invivo_D1_sig$NAME)
GSEA_invivo_D1_sig
nrow(GSEA_invivo_D1_sig)

# GSEA in vivo acute injury D3 vs Healthy
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/GSEA/GSEA_Results/20231208_invivo_D3.Gsea.1702042580607/')

GSEA_invivo_D3_up <- read.table('gsea_report_for_D3_1702042580607.tsv', header = TRUE, sep = '\t')
nrow(GSEA_invivo_D3_up)
GSEA_invivo_D3_up$Regulation <- rep(c('UP'), 31)
GSEA_invivo_D3_down <- read.table('gsea_report_for_D0_1702042580607.tsv', header = TRUE, sep = '\t')
nrow(GSEA_invivo_D3_down)
GSEA_invivo_D3_down$Regulation <- rep(c('DOWN'), 19)
GSEA_invivo_D3 <- rbind(GSEA_invivo_D3_up, GSEA_invivo_D3_down)
GSEA_invivo_D3 <- GSEA_invivo_D3[, colnames(GSEA_invivo_D3) %in% c('NAME', 'SIZE', 'ES', 'NES', 'NOM.p.val', 'FDR.q.val', 
                                                                   'Regulation')]
head(GSEA_invivo_D3)
tail(GSEA_invivo_D3)
GSEA_invivo_D3_sig <- GSEA_invivo_D3[GSEA_invivo_D3$FDR.q.val < 0.01, ]
GSEA_invivo_D3_sig$NAME <- gsub('HALLMARK_', '', GSEA_invivo_D3_sig$NAME)
GSEA_invivo_D3_sig
nrow(GSEA_invivo_D3_sig)

# GSEA in vivo acute injury D7 vs Healthy
setwd('D:/foldername/1. mCell_Dedifferentiation/RawData/GSEA/GSEA_Results/20231208_invivo_D7.Gsea.1702042605058/')

GSEA_invivo_D7_up <- read.table('gsea_report_for_D7_1702042605058.tsv', header = TRUE, sep = '\t')
nrow(GSEA_invivo_D7_up)
GSEA_invivo_D7_up$Regulation <- rep(c('UP'), 27)
GSEA_invivo_D7_down <- read.table('gsea_report_for_D0_1702042605058.tsv', header = TRUE, sep = '\t')
nrow(GSEA_invivo_D7_down)
GSEA_invivo_D7_down$Regulation <- rep(c('DOWN'), 23)
GSEA_invivo_D7 <- rbind(GSEA_invivo_D7_up, GSEA_invivo_D7_down)
GSEA_invivo_D7 <- GSEA_invivo_D7[, colnames(GSEA_invivo_D7) %in% c('NAME', 'SIZE', 'ES', 'NES', 'NOM.p.val', 'FDR.q.val', 
                                                                   'Regulation')]
head(GSEA_invivo_D7)
tail(GSEA_invivo_D7)
GSEA_invivo_D7_sig <- GSEA_invivo_D7[GSEA_invivo_D7$FDR.q.val < 0.01, ]
GSEA_invivo_D7_sig$NAME <- gsub('HALLMARK_', '', GSEA_invivo_D7_sig$NAME)
GSEA_invivo_D7_sig
nrow(GSEA_invivo_D7_sig)

# find the combinations
GSEA_invivo_D1 <- GSEA_invivo_D1[, 'NAME', 'FDR.q.val']





