library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomicRanges)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--bamfolder", help="the dir of bam")
argv <- add_argument(argv,"--win", help="the window length,default length 5000k",default=5000)
argv <- add_argument(argv,"--species", help="optional homo_sapiens and mus_musculus")
argv <- parse_args(argv)

#read args
bamfolder = argv$bamfolder
win = argv$win
species = argv$species

if (species == 'homo_sapiens'){
    hgref='hg38'
    } else if (species=='mus_musculus'){
    hgref='mm10'
}

bamFile <- list.files(bamfolder, pattern = '*.sorted.bam$')
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, resolution = win,hgref = hgref)

# add X,Y to ref_raw
ref_raw <- bambedObj$ref
options(scipen = 20)
length_X <- ref_raw@seqinfo@seqlengths[(length(ref_raw@seqinfo)-1)]
length_Y <- ref_raw@seqinfo@seqlengths[length(ref_raw@seqinfo)]
length_win <- win*1000

start_list <- seq(1,round(length_X/length_win)*length_win+1,length_win)
end_list <- c(seq(length_win,length_X,length_win),length_X)
ref_raw <- c(ref_raw,GRanges(seq="chrX",ranges = IRanges(start = start_list,end = end_list)))
start_list <- seq(1,round(length_Y/length_win)*length_win+1,length_win)
end_list <- c(seq(length_win,length_Y,length_win),length_Y)
ref_raw <- c(ref_raw,GRanges(seq="chrY",ranges = IRanges(start = start_list,end = end_list)))

bambedObj$ref <- ref_raw

gc <- get_gc(ref_raw,hgref = hgref)

# Getting raw read depth
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, seq = 'paired-end', hgref = hgref)
Y_raw <- coverageObj$Y
QCmetric_raw <- get_samp_QC(bambedObj)
ref_raw$gc <- gc

qcObj <- perform_qc(Y_raw = Y_raw, 
                    sampname_raw = sampname_raw, ref_raw = ref_raw, 
                    QCmetric_raw = QCmetric_raw)

if (length(qcObj$sampname) == 1){
    Y <- matrix(qcObj$Y,dimnames=list(paste0(qcObj$ref@seqnames,":",qcObj$ref@ranges),qcObj$sampname))
    } else {
    Y <- qcObj$Y
}
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric

# get gini coefficient for each cell
Gini <- get_gini(Y)

ref_raw <- ref_raw[which(paste0(ref_raw@seqnames,":",ref_raw@ranges) %in% row.names(Y))]

normal_gini <- max(Gini)

# first-pass CODEX2 run with no latent factors
normObj <- normalize_codex2_ns_noK(Y_qc = Y,gc_qc = ref_raw$gc,norm_index = which(Gini<=normal_gini))
# Ploidy initialization
ploidy <- initialize_ploidy(Y = Y, Yhat = normObj$Yhat, ref = ref_raw)
normObj.scope <- normalize_scope_foreach(Y_qc = Y, gc_qc = ref_raw$gc,
                                                K = 1, ploidyInt = ploidy,
                                                norm_index = which(Gini<=normal_gini), T = 1,
                                                beta0 = normObj$beta.hat, nCores = 2)

Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]
fGC.hat <- normObj.scope$fGC.hat[[which.max(normObj.scope$BIC)]]

chrs <- unique(as.character(seqnames(ref_raw)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs
for (chri in chrs) {
    message('\n', chri, '\n')
    segment_cs[[chri]] <- segment_CBScs(Y = Y,
                                    Yhat = Yhat,
                                    sampname = colnames(Y_raw),
                                    ref = ref_raw,
                                    chr = chri,
                                    mode = "integer", max.ns = 1)
}
iCN_sim <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))

write.table(iCN_sim,file="./cnv.tsv",sep="\t",row.names=FALSE)
write.table(data.frame(ref_raw@seqnames,ref_raw@ranges),file="./ref.tsv",sep="\t",row.names=FALSE)