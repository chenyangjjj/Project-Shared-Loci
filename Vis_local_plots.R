rm(list = ls())
library(dplyr)
library(ggplot2)
library(data.table)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(cowplot)
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


im = 1
SNP_loci = data_conjFDR$SNP[im]
CHR_loci = data_conjFDR$chrnum[im]
BP_floor = data_conjFDR$chrpos[im] - 250000
BP_ceiling = data_conjFDR$chrpos[im] + 250000
trait3.loc <- locus(data = data_trait3, ens_db = "EnsDb.Hsapiens.v75",
                    chrom = "CHR", pos = "BP", p = "PVAL",
                    index_snp = SNP_loci, seqname = CHR_loci,
                    xrange = c(BP_floor, BP_ceiling))
trait3.loc <- link_LD(trait3.loc, token = "8e7c11adaa33")
x = gg_scatter(trait3.loc,size = 4) +
    theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.line = element_line(size = 1, colour = "black")
    )
gg_addgenes(x, trait3.loc)
ggsave("/home/holstegelab-cjiang/projects/project_AD_shared/results/plots_rebuttal/gg_scatter.png",x)

ggsave("/home/holstegelab-cjiang/projects/project_AD_shared/results/plots_rebuttal/locus_ggplot,png",locus_ggplot(trait3.loc))
