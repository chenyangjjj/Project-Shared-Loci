rm(list = ls())
library(dplyr)
library(ggplot2)
library(data.table)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(cowplot)

main_path = "/home/holstegelab-cjiang/projects/project_AD_shared/"
# conjFDR results
file_conjFDR="/home/holstegelab-cjiang/projects/project_AD_shared/results/AD2023_loci.table"
data_conjFDR = fread(file_conjFDR)
Genome = fread("/home/holstegelab-cjiang/projects/admixture/ref_data/all_phase3.bim")
data_conjFDR = merge(data_conjFDR,Genome[,c("V1","V4","V2")],by.x = c("chrnum","chrpos"),by.y = c("V1","V4"),all.x = TRUE)


# cross-traits GWAS
data_trait1 <- fread("/home/holstegelab-cjiang/projects/project_AD_shared/data/AD_Bellenguz2023_DF_V1.csv")
data_trait2 <- fread("/home/holstegelab-cjiang/projects/project_AD_shared/data/conjFDR_Hipp_lh_biobank_clean.csv")
data_trait3 <- fread("/home/holstegelab-cjiang/projects/project_AD_shared/data/conjFDR_Hipp_rh_biobank_clean.csv")


for (im in c(6,2,10)) {

SNP_loci = data_conjFDR$SNP[im]
CHR_loci = data_conjFDR$chrnum[im]
BP_floor = data_conjFDR$chrpos[im] - 250000
BP_ceiling = data_conjFDR$chrpos[im] + 250000

trait1.loc <- locus(data = data_trait1, ens_db = "EnsDb.Hsapiens.v75",
                    chrom = "CHR", pos = "BP", p = "PVAL",
                    index_snp = SNP_loci, seqname = CHR_loci,
                    xrange = c(BP_floor, BP_ceiling))
trait1.loc <- link_LD(trait1.loc, token = "8e7c11adaa33")
locus_1 <- gg_scatter(trait1.loc,size = 5,pcutoff = 5e-20) +
    theme(
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 35),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 35),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 35),
    axis.line = element_line(size = 1, colour = "black")
    )
track1 = gg_genetracks(trait1.loc,cex.text = 0.7) +
            theme(
                legend.text = element_text(size = 35),
                legend.title = element_text(size = 35),
                axis.title.x = element_text(size = 35),
                axis.text.x = element_text(size = 35),
                axis.line = element_line(size = 1, colour = "black")
                )

trait2.loc <- locus(data = data_trait2, ens_db = "EnsDb.Hsapiens.v75",
                    chrom = "CHR", pos = "BP", p = "PVAL",
                    index_snp = SNP_loci, seqname = CHR_loci,
                    xrange = c(BP_floor, BP_ceiling))
trait2.loc <- link_LD(trait2.loc, token = "8e7c11adaa33")
locus_2 = gg_scatter(trait2.loc,size = 5,pcutoff = 5e-20) +
    theme(
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 35),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 35),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 35),
    axis.line = element_line(size = 1, colour = "black")
    )
track2 = gg_genetracks(trait2.loc,cex.text = 0.7) +
            theme(
                legend.text = element_text(size = 35),
                legend.title = element_text(size = 35),
                axis.title.x = element_text(size = 35),
                axis.text.x = element_text(size = 35),
                axis.line = element_line(size = 1, colour = "black")
                )

trait3.loc <- locus(data = data_trait3, ens_db = "EnsDb.Hsapiens.v75",
                    chrom = "CHR", pos = "BP", p = "PVAL",
                    index_snp = SNP_loci, seqname = CHR_loci,
                    xrange = c(BP_floor, BP_ceiling))
trait3.loc <- link_LD(trait3.loc, token = "8e7c11adaa33")
locus_3 = gg_scatter(trait3.loc,size = 5,pcutoff = 5e-20) +
    theme(
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 35),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 35),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 35),
    axis.line = element_line(size = 1, colour = "black")
    )
track3 = gg_genetracks(trait3.loc,cex.text = 0.8) +
            theme(
                legend.text = element_text(size = 35),
                legend.title = element_text(size = 35),
                axis.title.x = element_text(size = 35),
                axis.text.x = element_text(size = 35),
                axis.line = element_line(size = 1, colour = "black")
                )

locus_ALL = ggarrange(locus_1,locus_2,locus_3,nrow = 3,ncol = 1)
ggsave(paste0("/home/holstegelab-cjiang/projects/project_AD_shared/results/plots_rebuttal/gg_scatter_",SNP_loci,".png"),locus_ALL,width = 10,height = 27)
ggsave(paste0("/home/holstegelab-cjiang/projects/project_AD_shared/results/plots_rebuttal/gg_scatter_",SNP_loci,".pdf"),locus_ALL,width = 10,height = 27,dpi = 350)

ggsave(paste0("/home/holstegelab-cjiang/projects/project_AD_shared/results/plots_rebuttal/gg_track_",SNP_loci,".pdf"),track3,width = 10,height = 3)


}
