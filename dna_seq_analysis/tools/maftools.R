library(maftools)
library(stringr)
library(deconstructSigs)
library(BSgenome)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- add_argument(argv,"--species", help="Ensembl species name.")
argv <- add_argument(argv,"--maf_file", help="maf file")
argv <- parse_args(argv)

#read args
outdir = argv$outdir
species = argv$species
maf_file <- argv$maf_file

vep.laml = read.maf(maf = maf_file)

laml = vep.laml
laml@data=laml@data[!grepl('^MT-',laml@data$Hugo_Symbol),]
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
mut = laml@data[laml@data$t_alt_count >= 5 &
                  laml@data$t_vaf >= 0.05, c("Hugo_Symbol",
                                             "Chromosome",
                                             "Start_Position",
                                             "Tumor_Sample_Barcode",
                                             "t_vaf")]
mut = as.data.frame(mut)
write.csv(mut,file = paste0(outdir,'/08.annotated/mut_gene.txt'))

png(paste0(outdir,'/08.annotated/plotmafSummary_vep.png'),res = 150,width = 1080,height = 1500)
plotmafSummary(maf = vep.laml,
                 rmOutlier = TRUE,
                 showBarcodes = T,
                 addStat = 'median',
                 dashboard = TRUE,
                 titvRaw = FALSE)
dev.off()


if (species == 'homo_sapiens'){
    library(BSgenome.Hsapiens.NCBI.GRCh38)
    hgref=BSgenome.Hsapiens.NCBI.GRCh38
    chrom_list <- c(1:22,'X','Y','MT')
    } else if (species=='mus_musculus'){
    library(BSgenome.Mmusculus.UCSC.mm39)
    hgref=BSgenome.Mmusculus.UCSC.mm39
    chrom_list <- c(1:19,'X','Y','MT')
}

maf=read.table(maf_file,header = T,sep = '\t',quote = "")
sample.mut.ref <- data.frame(Sample=maf[,16], 
                            chr = maf[,5],
                            pos = maf[,6],
                            ref = maf[,11],
                            alt = maf[,13])
sample.mut.ref <- subset(sample.mut.ref, chr %in% chrom_list )

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = hgref)

df = data.frame()
for (i in rownames(sigs.input)){
    png(paste0(outdir,'/11.split/',i,'/plotSignature.png'),res = 150,width=1080,height=800)
    sigs.output <- whichSignatures(tumor.ref = sigs.input,
                                signatures.ref = signatures.cosmic, 
                                sample.id = i,
                                contexts.needed = TRUE)
    plotSignatures(sigs.output)
    dev.off()
    df = rbind(df,sigs.output$weights)
}

df = df[ , apply(df, 2, function(x){sum(x>0)})>0]
pheatmap::pheatmap(df,cluster_cols = F,cluster_rows = F,fontsize = 20,filename=paste0(outdir,'/08.annotated/PlotSignature_heapmap.png'))