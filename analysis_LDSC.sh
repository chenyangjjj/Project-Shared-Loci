
python munge_sumstats.py \
  --sumstats /home/holstegelab-cjiang/projects/project_AD_shared/data/AD_Bellenguz2023_DF_V1.csv \
  --out /home/holstegelab-cjiang/projects/project_AD_shared/rebuttal/data/AD_ldsc.sumstats \
  --ignore Z \
  --p PVAL \
  --merge-alleles w_hm3.snplist

python munge_sumstats.py \
  --sumstats /home/holstegelab-cjiang/projects/project_AD_shared/data/conjFDR_Hipp_lh_biobank_clean.csv \
  --out /home/holstegelab-cjiang/projects/project_AD_shared/rebuttal/data/Hipp_lh_ldsc.sumstats \
  --p PVAL \
  --signed-sumstats Z \
  --merge-alleles w_hm3.snplist

python munge_sumstats.py \
  --sumstats /home/holstegelab-cjiang/projects/project_AD_shared/data/conjFDR_Hipp_rh_biobank_clean.csv \
  --out /home/holstegelab-cjiang/projects/project_AD_shared/rebuttal/data/Hipp_rh_ldsc.sumstats \
  --p PVAL \
  --signed-sumstats Z \
  --merge-alleles w_hm3.snplist

  
python ldsc.py \
  --rg AD_ldsc.sumstats.gz,Hipp_rh_ldsc.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out output/genetic_correlation
