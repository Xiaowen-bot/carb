
#-------------------------------------------------------------------------------
#  December 5, 2024
#  The purpose of this program is to derive Carb's metabolites signature from LVS data
#  Programmer: Xiaowen Wang /udd/n2xwa/carb_met/
#-------------------------------------------------------------------------------

# devtools::install("/proj/nhmngs/nhmng00/R/chanmetab", build_vignettes = TRUE,
#                  upgrade = FALSE)
library(chanmetab)
library(Biobase)
# install.packages("stringr")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("MASS")
# install.packages("statar", lib = "/udd/n2xwa/R/x86_64-pc-linux-gnu-library/4.2", repos = "https://cran.r-project.org")
# install.packages("naniar")
# install.packages("ggpubr")
# install.packages("ggsci")
# install.packages("glmnet")
# install.packages("quickReg")
# install.packages("table1")

library(glmnet)
library(data.table)
library(cowplot)
# library(plotROC)
library(dplyr)
# library(quickReg)

library(ggrepel)
library(stringr)
library(ggplot2)
library(dplyr)
library(MASS)
library(haven)
.libPaths("/udd/n2xwa/R/x86_64-pc-linux-gnu-library/4.2")
library(statar)
library(table1)
library(plotROC)

setwd("/udd/n2xwa/carb_met")

all_lvs <- merge_metab_data(endpoints = c("lvs"),
                                collection_to_use = "substudy",
                                transformation = "transform_ln_z_score",
                                keep_unknown_metabolites = F,
                                keep_failed_pm_metabolites = F, 
                                impute_missing_function = "impute_one_half_min",
                                impute_cutoff=0.3,
                                combine_cohorts = TRUE)

table(all_lvs$summary$cohort) 
dim(all_lvs$expr_set$all_cohorts) # Features 405   Samples 1201


p_lvs=pData(all_lvs$expr_set$all_cohorts)
f_lvs=fData(all_lvs$expr_set$all_cohorts)
e_lvs=data.frame(exprs(all_lvs$expr_set$all_cohorts))

exprs(all_lvs$expr_set$all_cohorts)
all_lvs$expr_set$all_cohorts


f_lvs_cv <- f_lvs %>% filter(mean_cv<30 & mean_icc>0.4) # 405 to 341 metabolites
e_lvs_cv <- subset(e_lvs, (rownames(e_lvs) %in% rownames(f_lvs_cv))) # 341 metabolites

a_lvs=data.frame(t(e_lvs_cv))
rownames(p_lvs) <- paste("X", rownames(p_lvs), sep="")
identical(rownames(a_lvs), rownames(p_lvs)) # should be TRUE
z_lvs=cbind(a_lvs,p_lvs) 
#head(z_lvs$id)
z_lvs=z_lvs %>%  mutate(id=str_sub(id, start = 1, end = 6),
                        newid=case_when(cohort=="nhs1" ~ as.numeric(id)+1000000,
                                        cohort=="nhs2" ~ as.numeric(id)+2000000,
                                        cohort=="hpfs" ~ as.numeric(id)+3000000))

non_na_counts <- colMeans(is.na(z_lvs[,1:341]))
z_lvs_part1 <- z_lvs[,1:341][, non_na_counts <0.3] # 341 to 295 metabolites
z_lvs=cbind(z_lvs_part1, z_lvs[, 342:364])
# library(naniar)
# missing <- miss_case_summary(z_lvs[,1:295])
# View(missing) # case 725 and 28 are still missing
z_lvs <-  z_lvs[-c(725,28),]  #1201 obs to 1199 obs

#combine AHEI data#

nhs_ahei <- read.csv("/udd/n2xwa/lvs_met/nhs_ahei.csv")
nhs2_ahei <- read.csv("/udd/n2xwa/lvs_met/nhs2_ahei.csv")
hpfs_ahei <- read.csv("/udd/n2xwa/lvs_met/hpfs_ahei.csv")

nhs_ahei <- nhs_ahei %>% 
  mutate(ahei_a = rowMeans(cbind(ahei86_a,ahei90_a,
                                     ahei94_a,ahei98_a,
                                     ahei02_a,ahei06_a,ahei10_a),na.rm = TRUE),
         ahei_nowgr = rowMeans(cbind(ahei86_nowgr,ahei90_nowgr,
                                     ahei94_nowgr,ahei98_nowgr,
                                     ahei02_nowgr,ahei06_nowgr,ahei10_nowgr), na.rm = TRUE),
         ahei_noveg = rowMeans(cbind(ahei86_noveg,ahei90_noveg,
                                     ahei94_noveg,ahei98_noveg,
                                     ahei02_noveg,ahei06_noveg,ahei10_noveg), na.rm = TRUE),
         ahei_nofrt = rowMeans(cbind(ahei86_nofrt,ahei90_nofrt,
                                     ahei94_nofrt,ahei98_nofrt,
                                     ahei02_nofrt,ahei06_nofrt,ahei10_nofrt), na.rm = TRUE),
         ahei_noleg = rowMeans(cbind(ahei86_noleg,ahei90_noleg,
                                     ahei94_noleg,ahei98_noleg,
                                     ahei02_noleg,ahei06_noleg,ahei10_noleg), na.rm = TRUE),
         ahei_na = rowMeans (cbind(ahei86_na,ahei90_na,
                                 ahei94_na,ahei98_na,
                                 ahei02_na,ahei06_na,ahei10_na), na.rm = TRUE)) %>%
  dplyr::select(id, ahei_a, ahei_nowgr, ahei_noveg,ahei_nofrt,ahei_noleg,ahei_na)

nhs2_ahei <- nhs2_ahei %>% 
  mutate(ahei_a = rowMeans(cbind(ahei91_a,ahei95_a,
                                  ahei99_a,ahei03_a,
                                  ahei07_a,ahei11_a), na.rm = TRUE),
         ahei_nowgr = rowMeans(cbind(ahei91_nowgr,ahei95_nowgr,
                                     ahei99_nowgr,ahei03_nowgr,
                                     ahei07_nowgr,ahei11_nowgr), na.rm = TRUE),
         ahei_noveg = rowMeans(cbind(ahei91_noveg,ahei95_noveg,
                                     ahei99_noveg,ahei03_noveg,
                                     ahei07_noveg,ahei11_noveg), na.rm = TRUE),
         ahei_nofrt = rowMeans(cbind(ahei91_nofrt,ahei95_nofrt,
                                     ahei99_nofrt,ahei03_nofrt,
                                     ahei07_nofrt,ahei11_nofrt), na.rm = TRUE),
         ahei_noleg = rowMeans(cbind(ahei91_noleg,ahei95_noleg,
                                     ahei99_noleg,ahei03_noleg,
                                     ahei07_noleg,ahei11_noleg), na.rm = TRUE),
         ahei_na = rowMeans(cbind(ahei91_na,ahei95_na,
                                 ahei99_na,ahei03_na,
                                 ahei07_na,ahei11_na), na.rm = TRUE)) %>%
  dplyr::select(id, ahei_a, ahei_nowgr,  ahei_noveg,ahei_nofrt,ahei_noleg, ahei_na)


hpfs_ahei <- hpfs_ahei %>%mutate(ahei_a = rowMeans(cbind(ahei86_a,ahei90_a,
                                                          ahei94_a,ahei98_a,
                                                          ahei02_a,ahei06_a,
                                                          ahei10_a), na.rm = TRUE),
                                 ahei_nowgr = rowMeans(cbind(ahei86_nowgr,ahei90_nowgr,
                                                             ahei94_nowgr,ahei98_nowgr,
                                                             ahei02_nowgr,ahei06_nowgr
                                                             ,ahei10_nowgr), na.rm = TRUE),
                                 ahei_noveg = rowMeans(cbind(ahei86_noveg,ahei90_noveg,
                                                             ahei94_noveg,ahei98_noveg,
                                                             ahei02_noveg,ahei06_noveg,ahei10_noveg), na.rm = TRUE),
                                 ahei_nofrt = rowMeans(cbind(ahei86_nofrt,ahei90_nofrt,
                                                             ahei94_nofrt,ahei98_nofrt,
                                                             ahei02_nofrt,ahei06_nofrt,ahei10_nofrt), na.rm = TRUE),
                                 ahei_noleg = rowMeans(cbind(ahei86_noleg,ahei90_noleg,
                                                             ahei94_noleg,ahei98_noleg,
                                                             ahei02_noleg,ahei06_noleg,ahei10_noleg), na.rm = TRUE),
                                 ahei_na = rowMeans(cbind(ahei86_na,ahei90_na,
                                                         ahei94_na,ahei98_na,
                                                         ahei02_na,ahei06_na,
                                                         ahei10_na), na.rm = TRUE)) %>%
  dplyr::select(id, ahei_a, ahei_nowgr, ahei_noveg,ahei_nofrt,ahei_noleg, ahei_na)

ahei<- rbind(nhs_ahei,nhs2_ahei,hpfs_ahei)


# Import MLVS data ####

# 7DDR food data
setwd("/proj/polvss/polvs00/MLVS/data_for_analysis/diet_7ddr/sas_data/")
data_mlvs_food <- haven::read_sas("dr_food_wkt_mean.sas7bdat")

data_mlvs_food <- data_mlvs_food %>% mutate(whole_grain=rowSums(cbind(cer07_wf_dr_wt, oat07_wf_dr_wt, ckcer07_wf_dr_wt, 
                                                                      ryebr07_wf_dr_wt, dkbr07_wf_dr_wt,brice07_wf_dr_wt, oatbr07_wf_dr_wt, 
                                                                      bran07_wf_dr_wt, ffpop07_wf_dr_wt, popc07_wf_dr_wt),na.rm = T),
                                            refined_grain=rowSums(cbind(whbr07_wf_dr_wt, crack07_wf_dr_wt,pizza07_wf_dr_wt,
                                                                        engl07_wf_dr_wt, muff07_wf_dr_wt, pcake07_wf_dr_wt,  
                                                                        wrice07_wf_dr_wt, pasta07_wf_dr_wt,pretz07_wf_dr_wt),na.rm = T),
                                            fruit=rowSums(cbind(rais07_wf_dr_wt,prune07_wf_dr_wt,prunj07_wf_dr_wt,
                                                                 ban07_wf_dr_wt,cant07_wf_dr_wt,avo07_wf_dr_wt,
                                                                 appl07_wf_dr_wt,aj07_wf_dr_wt,oran07_wf_dr_wt,
                                                                 ojca07_wf_dr_wt,oj07_wf_dr_wt,grfr07_wf_dr_wt,
                                                                 othj07_wf_dr_wt,straw07_wf_dr_wt,blueb07_wf_dr_wt,
                                                                 peach07_wf_dr_wt,apric07_wf_dr_wt),na.rm = T),
                                            vegetable=rowSums(cbind(tom07_wf_dr_wt,rtoj0ywkt,rtosau0ywkt,salsa07_wf_dr_wt,rsbean07wkt,
                                                                    bean07_wf_dr_wt,tofu07_wf_dr_wt,peas07_wf_dr_wt,brocc07_wf_dr_wt,
                                                                    cauli07_wf_dr_wt,cabb07_wf_dr_wt,bruss07_wf_dr_wt,rcar07_wf_dr_wt,
                                                                    ccar07_wf_dr_wt,corn07_wf_dr_wt,mixv07_wf_dr_wt,yam07_wf_dr_wt,
                                                                    osqua07_wf_dr_wt,eggpl07_wf_dr_wt,kale07_wf_dr_wt,cspin07_wf_dr_wt,
                                                                    rspin07_wf_dr_wt,ilett07_wf_dr_wt,rlett07_wf_dr_wt,cel07_wf_dr_wt,
                                                                    grpep07_wf_dr_wt,oniog07_wf_dr_wt,oniov07_wf_dr_wt),na.rm = T),
                                            legume=rowSums(cbind(rsbean07wkt,bean07_wf_dr_wt,tofu07_wf_dr_wt,peas07_wf_dr_wt),na.rm = T),
                                            SSB=rowSums(cbind(cola07_wf_dr_wt, otsug07_wf_dr_wt, punch07_wf_dr_wt),na.rm = T),
                                            nonstarchyveg=rowSums(cbind(tom07_wf_dr_wt,rtoj0ywkt,rtosau0ywkt,salsa07_wf_dr_wt,rsbean07wkt,
                                                                        bean07_wf_dr_wt,tofu07_wf_dr_wt,brocc07_wf_dr_wt,
                                                                        cauli07_wf_dr_wt,cabb07_wf_dr_wt,bruss07_wf_dr_wt,rcar07_wf_dr_wt,
                                                                        ccar07_wf_dr_wt,mixv07_wf_dr_wt,
                                                                        osqua07_wf_dr_wt,eggpl07_wf_dr_wt,kale07_wf_dr_wt,cspin07_wf_dr_wt,
                                                                        rspin07_wf_dr_wt,ilett07_wf_dr_wt,rlett07_wf_dr_wt,cel07_wf_dr_wt,
                                                                        grpep07_wf_dr_wt,oniog07_wf_dr_wt,oniov07_wf_dr_wt),na.rm = T),
                                            starchyveg=rowSums(cbind(fries07_wf_dr_wt,pot07_wf_dr_wt,pchip07_wf_dr_wt,
                                                                     yam07_wf_dr_wt,peas07_wf_dr_wt,corn07_wf_dr_wt),na.rm = T),    
                                            potato=rowSums(cbind(fries07_wf_dr_wt,pot07_wf_dr_wt,pchip07_wf_dr_wt),na.rm = T),                                                                                                   
                                            cohort='mlvs') %>% dplyr::select(id, whole_grain, refined_grain,fruit,vegetable,legume,nonstarchyveg,starchyveg,potato,SSB,cohort)
# summary(data_mlvs_food)

# 7ddr nutrients data 

setwd("/proj/polvss/polvs00/MLVS/data_for_analysis/diet_7ddr/sas_data/")
data_mlvs_nts <- haven::read_sas("dr_nts_wkt_mean.sas7bdat")
data_mlvs_nts<-data_mlvs_nts %>% mutate(cohort="mlvs") %>% dplyr::select(id, pctcho_fo_dr_wtavg,aofib_fo_dr_wtavg,
                                                                         pctfat_fo_dr_wtavg,pctpro_fo_dr_wtavg,calor_fo_dr_wtavg,
                                                                         addsug_fo_dr_wtavg,cohort,asugar_fo_dr_wtavg,gluten_fo_dr_wtavg)
# summary(data_mlvs_nts)


# Import WLVS data ####


# 7DDR food data

setwd("/proj/polvss/polvs00/WLVS/data_for_analysis/diet_7ddr/sas_data/")
data_wlvs_food <- haven::read_sas("dr_food_wkt_mean.sas7bdat")
data_wlvs_food <- data_wlvs_food %>% mutate(whole_grain=rowSums(cbind(cer07_wf_dr_wt, oat07_wf_dr_wt, ckcer07_wf_dr_wt, 
                                                                      ryebr07_wf_dr_wt, dkbr07_wf_dr_wt,brice07_wf_dr_wt, oatbr07_wf_dr_wt, 
                                                                      bran07_wf_dr_wt, ffpop07_wf_dr_wt, popc07_wf_dr_wt),na.rm = T),
                                            refined_grain=rowSums(cbind(whbr07_wf_dr_wt, crack07_wf_dr_wt, 
                                                                      engl07_wf_dr_wt, muff07_wf_dr_wt, pcake07_wf_dr_wt,pizza07_wf_dr_wt,
                                                                      wrice07_wf_dr_wt, pasta07_wf_dr_wt,pretz07_wf_dr_wt),na.rm = T),
                                            fruit=rowSums(cbind(rais07_wf_dr_wt,prune07_wf_dr_wt,prunj07_wf_dr_wt,
                                                              ban07_wf_dr_wt,cant07_wf_dr_wt,avo07_wf_dr_wt,
                                                              appl07_wf_dr_wt,aj07_wf_dr_wt,oran07_wf_dr_wt,
                                                              ojca07_wf_dr_wt,oj07_wf_dr_wt,grfr07_wf_dr_wt,
                                                              othj07_wf_dr_wt,straw07_wf_dr_wt,blueb07_wf_dr_wt,
                                                              peach07_wf_dr_wt,apric07_wf_dr_wt),na.rm = T),
                                            vegetable=rowSums(cbind(tom07_wf_dr_wt,rtoj0ywkt,rtosau0ywkt,salsa07_wf_dr_wt,rsbean07wkt,
                                                                  bean07_wf_dr_wt,tofu07_wf_dr_wt,peas07_wf_dr_wt,brocc07_wf_dr_wt,
                                                                  cauli07_wf_dr_wt,cabb07_wf_dr_wt,bruss07_wf_dr_wt,rcar07_wf_dr_wt,
                                                                  ccar07_wf_dr_wt,corn07_wf_dr_wt,mixv07_wf_dr_wt,yam07_wf_dr_wt,
                                                                  osqua07_wf_dr_wt,eggpl07_wf_dr_wt,kale07_wf_dr_wt,cspin07_wf_dr_wt,
                                                                  rspin07_wf_dr_wt,ilett07_wf_dr_wt,rlett07_wf_dr_wt,cel07_wf_dr_wt,
                                                                  grpep07_wf_dr_wt,oniog07_wf_dr_wt,oniov07_wf_dr_wt),na.rm = T),
                                           legume=rowSums(cbind(rsbean07wkt,bean07_wf_dr_wt,tofu07_wf_dr_wt,peas07_wf_dr_wt),na.rm = T),
                                           SSB=rowSums(cbind(cola07_wf_dr_wt, otsug07_wf_dr_wt, punch07_wf_dr_wt),na.rm = T),
                                           nonstarchyveg=rowSums(cbind(tom07_wf_dr_wt,rtoj0ywkt,rtosau0ywkt,salsa07_wf_dr_wt,rsbean07wkt,
                                                                      bean07_wf_dr_wt,tofu07_wf_dr_wt,brocc07_wf_dr_wt,
                                                                      cauli07_wf_dr_wt,cabb07_wf_dr_wt,bruss07_wf_dr_wt,rcar07_wf_dr_wt,
                                                                      ccar07_wf_dr_wt,mixv07_wf_dr_wt,
                                                                      osqua07_wf_dr_wt,eggpl07_wf_dr_wt,kale07_wf_dr_wt,cspin07_wf_dr_wt,
                                                                      rspin07_wf_dr_wt,ilett07_wf_dr_wt,rlett07_wf_dr_wt,cel07_wf_dr_wt,
                                                                      grpep07_wf_dr_wt,oniog07_wf_dr_wt,oniov07_wf_dr_wt),na.rm = T),
                                           starchyveg=rowSums(cbind(fries07_wf_dr_wt,pot07_wf_dr_wt,pchip07_wf_dr_wt,
                                                                   yam07_wf_dr_wt,peas07_wf_dr_wt,corn07_wf_dr_wt),na.rm = T),
                                           potato=rowSums(cbind(fries07_wf_dr_wt,pot07_wf_dr_wt,pchip07_wf_dr_wt),na.rm = T), 
                                           cohort='wlvs') %>% dplyr::select(id, whole_grain, refined_grain,fruit,vegetable,legume,nonstarchyveg,starchyveg,potato,SSB,cohort)
# summary(data_wlvs_food)


# 7ddr nutrients data #

setwd("/proj/polvss/polvs00/WLVS/data_for_analysis/diet_7ddr/sas_data/")
data_wlvs_nts <- haven::read_sas("dr_nts_wkt_mean.sas7bdat")
data_wlvs_nts <- data_wlvs_nts %>% mutate(cohort="wlvs") %>% dplyr::select(id, pctcho_fo_dr_wtavg, aofib_fo_dr_wtavg,
                                                                           pctfat_fo_dr_wtavg,pctpro_fo_dr_wtavg,
                                                                           calor_fo_dr_wtavg,addsug_fo_dr_wtavg,cohort,asugar_fo_dr_wtavg,gluten_fo_dr_wtavg)

# summary(data_wlvs_nts)


## combine LVS data ####

data_lvs_food<- rbind(data_mlvs_food,data_wlvs_food)# 1466
data_lvs_nts <- rbind(data_mlvs_nts, data_wlvs_nts)# 1466
table(data_lvs_food$cohort)
table(data_lvs_nts$cohort)
# mlvs wlvs 
# 692  774 


# MLVS ####

# FFQ data ####
setwd("/proj/polvss/polvs00/MLVS/data_for_analysis/diet_ffq/sas_data/")
data_mlvs_ffq1_nts <- haven::read_sas("mlvsffq1_excl_nts.sas7bdat") 

data_mlvs_ffq1_nts<-data_mlvs_ffq1_nts  %>% dplyr::select(id,calor_fs_ffq1,carbo_fs_ffq1,frtcb_njs_fs_ffq1,
                                                          leg_carb_fs_ffq1,pot_carb_fs_ffq1,ref_grn_fs_ffq1,
                                                          sum_whgcarb_fs_ffq1,vegcb_npl_fs_ffq1,
                                                          addsug_fs_ffq1,adds_fds_fs_ffq1)

data_mlvs_ffq2_nts <- haven::read_sas("mlvsffq2_excl_nts.sas7bdat") 

data_mlvs_ffq2_nts<-data_mlvs_ffq2_nts  %>% dplyr::select(id,calor_fs_ffq2,carbo_fs_ffq2,frtcb_njs_fs_ffq2,
                                                          leg_carb_fs_ffq2,pot_carb_fs_ffq2,ref_grn_fs_ffq2,
                                                          sum_whgcarb_fs_ffq2,vegcb_npl_fs_ffq2,
                                                          addsug_fs_ffq2,adds_fds_fs_ffq2)
combined_mlvs<-merge(data_mlvs_ffq1_nts, data_mlvs_ffq2_nts, by = 'id',all = FALSE, sort = TRUE)


ffq_mlvs <- combined_mlvs %>%
  rowwise() %>%
  mutate(cohort = "mlvs",
         avg_calor_fs = mean(c(calor_fs_ffq1, calor_fs_ffq2), na.rm = TRUE),
         avg_carbo_fs = mean(c(carbo_fs_ffq1, carbo_fs_ffq2), na.rm = TRUE),
         avg_frtcb_njs_fs = mean(c(frtcb_njs_fs_ffq1, frtcb_njs_fs_ffq2), na.rm = TRUE),
         avg_leg_carb_fs = mean(c(leg_carb_fs_ffq1, leg_carb_fs_ffq2), na.rm = TRUE),
         avg_pot_carb_fs = mean(c(pot_carb_fs_ffq1, pot_carb_fs_ffq2), na.rm = TRUE),
         avg_ref_grn_fs = mean(c(ref_grn_fs_ffq1, ref_grn_fs_ffq2), na.rm = TRUE),
         avg_sum_whgcarb_fs = mean(c(sum_whgcarb_fs_ffq1, sum_whgcarb_fs_ffq2), na.rm = TRUE),
         avg_vegcb_npl_fs = mean(c(vegcb_npl_fs_ffq1, vegcb_npl_fs_ffq2), na.rm = TRUE),
         avg_addsug_fs = mean(c(addsug_fs_ffq1, addsug_fs_ffq2), na.rm = TRUE),
         avg_adds_fds_fs = mean(c(adds_fds_fs_ffq1, adds_fds_fs_ffq2), na.rm = TRUE)
  ) %>% ungroup() %>% dplyr::select(id, cohort,starts_with("avg_"))


# WLVS ####

# FFQ data ####
setwd("/proj/polvss/polvs00/WLVS/data_for_analysis/diet_ffq/sas_data/")
data_wlvs_ffq1_nts <- haven::read_sas("wlvsffq1_excl_nts.sas7bdat") 

data_wlvs_ffq1_nts<-data_wlvs_ffq1_nts  %>% dplyr::select(id,calor_fs_ffq1,carbo_fs_ffq1,frtcb_njs_fs_ffq1,
                                                          leg_carb_fs_ffq1,pot_carb_fs_ffq1,ref_grn_fs_ffq1,
                                                          sum_whgcarb_fs_ffq1,vegcb_npl_fs_ffq1,
                                                          addsug_fs_ffq1,adds_fds_fs_ffq1)

data_wlvs_ffq2_nts <- haven::read_sas("wlvsffq2_excl_nts.sas7bdat") 

data_wlvs_ffq2_nts<-data_wlvs_ffq2_nts  %>% dplyr::select(id,calor_fs_ffq2,carbo_fs_ffq2,frtcb_njs_fs_ffq2,
                                                          leg_carb_fs_ffq2,pot_carb_fs_ffq2,ref_grn_fs_ffq2,
                                                          sum_whgcarb_fs_ffq2,vegcb_npl_fs_ffq2,
                                                          addsug_fs_ffq2,adds_fds_fs_ffq2)
combined_wlvs<-merge(data_wlvs_ffq1_nts, data_wlvs_ffq2_nts, by = 'id',all = FALSE, sort = TRUE)


ffq_wlvs <- combined_wlvs %>%
  rowwise() %>%
  mutate(cohort = "wlvs",
         avg_calor_fs = mean(c(calor_fs_ffq1, calor_fs_ffq2), na.rm = TRUE),
         avg_carbo_fs = mean(c(carbo_fs_ffq1, carbo_fs_ffq2), na.rm = TRUE),
         avg_frtcb_njs_fs = mean(c(frtcb_njs_fs_ffq1, frtcb_njs_fs_ffq2), na.rm = TRUE),
         avg_leg_carb_fs = mean(c(leg_carb_fs_ffq1, leg_carb_fs_ffq2), na.rm = TRUE),
         avg_pot_carb_fs = mean(c(pot_carb_fs_ffq1, pot_carb_fs_ffq2), na.rm = TRUE),
         avg_ref_grn_fs = mean(c(ref_grn_fs_ffq1, ref_grn_fs_ffq2), na.rm = TRUE),
         avg_sum_whgcarb_fs = mean(c(sum_whgcarb_fs_ffq1, sum_whgcarb_fs_ffq2), na.rm = TRUE),
         avg_vegcb_npl_fs = mean(c(vegcb_npl_fs_ffq1, vegcb_npl_fs_ffq2), na.rm = TRUE),
         avg_addsug_fs = mean(c(addsug_fs_ffq1, addsug_fs_ffq2), na.rm = TRUE),
         avg_adds_fds_fs = mean(c(adds_fds_fs_ffq1, adds_fds_fs_ffq2), na.rm = TRUE)
  ) %>% ungroup() %>% dplyr::select(id, cohort,starts_with("avg_"))

# LVS FFQ data
ffq_lvs<- rbind(ffq_mlvs,ffq_wlvs)# 1367


## LVS combine: combine metabolites, adjusted-ahei, ffq_lvs ####

temp_lvs_t<-merge(z_lvs, data_lvs_nts[,-8], by = 'id',all = FALSE, sort = TRUE) # 1196 obs. 
ahei$id <- as.character(ahei$id)
temp_lvs_t$id <- as.character(temp_lvs_t$id)
temp_lvs_t <- left_join(temp_lvs_t, ahei, by = "id") 
temp_lvs_t<-temp_lvs_t %>% filter(!duplicated(id)) # 1196 obs.

ffq_lvs$id <- as.character(ffq_lvs$id)
temp_lvs <- left_join(temp_lvs_t, ffq_lvs[ , -2], by = 'id')
temp_lvs<-temp_lvs %>% filter(!duplicated(id)) # 1196 obs.



##  ####

temp_lvs<-temp_lvs %>% group_by(cohort) %>% mutate(white=case_when(race=='White' ~ 1,T ~ 0),
                                                   fasting=case_when(fast=='nonfasting' ~ 0,T ~ 1),
                                                   smkstatus=case_when(smoke=='current'~ 'current',
                                                                       smoke=='past' ~ 'past',
                                                                       T ~ 'never'),
                                                   neversmoking=case_when(smkstatus!='never' ~ 'no',
                                                                          T ~ 'yes'),
                                                   sex=case_when(cohort=='hpfs' ~ 'male',
                                                                 T ~ 'female'),
                                                   bmi_group=case_when(bmi<25 ~ 'normal',
                                                                       bmi<30 ~ 'overweight',
                                                                       bmi>=30 ~ 'obese',
                                                                       T ~ 'NA'),
                                                   act_bld = tidyr::replace_na(act, median(act, na.rm = TRUE)),
                                                   alco_bld = tidyr::replace_na(alco, median(alco, na.rm = TRUE)),
                                                   ahei_av = tidyr::replace_na(ahei_a, median(ahei_a, na.rm = TRUE)),
                                                   ahei_av_nowgr = tidyr::replace_na(ahei_nowgr, median(ahei_nowgr, na.rm = TRUE)),
                                                   ahei_av_noveg = tidyr::replace_na(ahei_noveg, median(ahei_noveg, na.rm = TRUE)),
                                                   ahei_av_nofrt = tidyr::replace_na(ahei_nofrt, median(ahei_nofrt, na.rm = TRUE)),
                                                   ahei_av_noleg = tidyr::replace_na(ahei_noleg, median(ahei_noleg, na.rm = TRUE)),
                                                   calories = tidyr::replace_na(calor_fo_dr_wtavg, median(calor_fo_dr_wtavg, na.rm = TRUE)),
                                                   per_carb = tidyr::replace_na(pctcho_fo_dr_wtavg, median(pctcho_fo_dr_wtavg, na.rm = TRUE)),
                                                   per_fat =  tidyr::replace_na(pctfat_fo_dr_wtavg, median(pctfat_fo_dr_wtavg, na.rm = TRUE)),
                                                   per_protein = tidyr::replace_na(pctpro_fo_dr_wtavg, median(pctpro_fo_dr_wtavg, na.rm = TRUE)),
                                                   addsugar = tidyr::replace_na(addsug_fo_dr_wtavg, median(addsug_fo_dr_wtavg, na.rm = TRUE)),
                                                   fiber = tidyr::replace_na(aofib_fo_dr_wtavg, median(aofib_fo_dr_wtavg, na.rm = TRUE)),
                                                   alco_g=statar::xtile(alco_bld, n = 5),
                                                   act_g=statar::xtile(act_bld, n = 5),
                                                   calor_g=statar::xtile(calories, n = 5),
                                                   prot_g=statar::xtile(per_protein, n = 5),
                                                   fat_g=statar::xtile(per_fat, n = 5),
                                                   ahei_g=statar::xtile(ahei_av, n = 5),
                                                   ahei_nowgr_g=statar::xtile(ahei_av_nowgr, n = 5),
                                                   ahei_noveg_g=statar::xtile(ahei_av_noveg, n = 5),
                                                   ahei_nofrt_g=statar::xtile(ahei_av_nofrt, n = 5),
                                                   ahei_noleg_g=statar::xtile(ahei_av_noleg, n = 5)) %>% as.data.frame()
names(data_lvs_food)

temp_lvs<-merge(temp_lvs, data_lvs_food[,-11], by = 'id',all = FALSE, sort = TRUE) # 1196

# derive energy-residual variables
model1 <- lm(refined_grain ~ calories, data = temp_lvs)
model2 <- lm(whole_grain ~ calories, data = temp_lvs)
model3 <- lm(fruit ~ calories, data = temp_lvs)
model4 <- lm(nonstarchyveg ~ calories, data = temp_lvs)
model5 <- lm(starchyveg ~ calories, data = temp_lvs)
model6 <- lm(addsugar ~ calories, data = temp_lvs)
model7 <- lm(vegetable ~ calories, data = temp_lvs)
model8 <- lm(legume ~ calories, data = temp_lvs)
model9 <- lm(potato ~ calories, data = temp_lvs)
model10 <- lm(gluten_fo_dr_wtavg ~ calories, data = temp_lvs)
# model8 <- lm(per_carb ~ fiber, data = temp_lvs)
# model9 <- lm(per_carb ~ SSB + refined_grain + potato, data = temp_lvs)

residuals1 <- resid(model1)
residuals2 <- resid(model2)
residuals3 <- resid(model3)

residuals4 <- resid(model4)
residuals5 <- resid(model5)
residuals6 <- resid(model6)
residuals7 <- resid(model7)
residuals8 <- resid(model8)
residuals9 <- resid(model9)
residuals10 <- resid(model10)

temp_lvs <- temp_lvs %>%
  mutate(refined_grain_r = residuals1,whole_grain_r=residuals2,fruit_r=residuals3,
         nonstarchyveg_r=residuals4,starchyveg_r=residuals5,addsugar_r=residuals6,
         vegetable_r=residuals7,legume_r=residuals8,potato_r=residuals9,gluten_r=residuals10)

# delete drug high-related mets: alpha-hydroxymetoprolol and acetaminophen
temp_lvs <- temp_lvs %>% 
  dplyr::select(-c(HMDB0060994, HMDB0001859))

names(temp_lvs)

save(temp_lvs, file = "/udd/n2xwa/carb_met/temp_lvs.RData")


#table 1#

T1 <- table1(~ageyr+sex+race+fasting+neversmoking+act_bld+alco_bld+
               bmi+bmi_group+calories+ahei_av+
               whole_grain+refined_grain+fruit+vegetable+starchyveg+nonstarchyveg+
               potato+legume+per_protein+per_carb++addsugar | cohort, data=temp_lvs, overall="Total")

print(T1)
write.csv(T1, "/udd/n2xwa/carb_met/table1.csv")

#-------------------------------------------------------------------------------

#                      Individual regression: metabolites~diet 

#-------------------------------------------------------------------------------

names(temp_lvs)
met_list=as.character(colnames(temp_lvs[,2:294]))
results_lvs=data.frame(met=met_list)
met_aligned=results_lvs %>% pull(met)

for (i in 1:293) {
  fit=lm(as.formula(paste0(met_aligned[i] ,"~ whole_grain_r + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_wgrain"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_wgrain"]=summ$coefficients[2,2]
  results_lvs[i, "p_wgrain"]=summ$coefficients[2,4]
  results_lvs[i, "n_wgrain"]=nrow(model.frame(fit))
  fit=lm(as.formula(paste0(met_aligned[i], "~ whole_grain_r + cohort + ageyr + white + factor(bmi_group) + factor(neversmoking) + 
  factor(ahei_nowgr_g) + factor(calor_g) + factor(act_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_wgrain_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_wgrain_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_wgrain_adj"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i],"~ addsugar_r + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_addsugar"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_addsugar"]=summ$coefficients[2,2]
  results_lvs[i, "p_addsugar"]=summ$coefficients[2,4]
  results_lvs[i, "n_addsugar"]=nrow(model.frame(fit))
  fit=lm(as.formula(paste0(met_aligned[i], "~ addsugar_r + cohort + ageyr  + white + factor(bmi_group) + factor(neversmoking) +
  factor(alco_g)+ factor(ahei_g) + factor(calor_g) + factor(act_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_addsugar_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_addsugar_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_addsugar_adj"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i],"~ per_carb + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_carb"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_carb"]=summ$coefficients[2,2]
  results_lvs[i, "p_carb"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ per_carb + cohort  + factor(bmi_group) + ageyr  + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                           factor(act_g) + factor(alco_g)+factor(prot_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_carb_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_carb_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_carb_adj"]=summ$coefficients[2,4]
  
}

results_lvs<- results_lvs %>% mutate(
  wgrain_intake=case_when(beta_wgrain_adj>0 & p_wgrain_adj<0.05/293  ~ 'positive',
                          beta_wgrain_adj<0 & p_wgrain_adj<0.05/293  ~ 'negative',
                          T ~ 'insignificant'),
  addsugar_intake=case_when(beta_addsugar_adj>0 & p_addsugar_adj<0.05/293  ~ 'positive',
                          beta_addsugar_adj<0 & p_addsugar_adj<0.05/293  ~ 'negative',
                          T ~ 'insignificant'),
  carbo_intake=case_when(beta_carb_adj>0 & p_carb_adj<0.05/293  ~ 'positive',
                         beta_carb_adj<0 & p_carb_adj<0.05/293  ~ 'negative',
                         T ~ 'insignificant'),
  wgrain_sig=case_when(p_wgrain_adj<0.05/293  ~ 100, T ~ 0),
  addsugar_sig=case_when(p_addsugar_adj<0.05/293  ~ 1, T ~ 0),
  carb_sig = case_when(p_carb_adj<0.05/293  ~ 10, T ~ 0),
  pattern1=wgrain_sig+addsugar_sig+carb_sig,
  pattern=case_when(pattern1==1 ~ 'addsugar', 
                    pattern1==10 ~ 'carb', 
                    pattern1==100 ~ 'wgrain', 
                    pattern1==11 ~ 'addsugar + carb', 
                    pattern1==101 ~ 'wgrain + addsugar', 
                    pattern1==110 ~ 'wgrain + carb',
                    pattern1==111 ~ 'wgrain + addsugar + carb', 
                    T ~ 'insignificant')) %>% dplyr::select(-c(wgrain_sig,addsugar_sig, carb_sig,pattern1))

table(results_lvs$pattern)


f_lvs_match<-f_lvs %>% mutate(met=rownames(f_lvs)) %>% 
  dplyr::select("met",'method',"metabolite_name",'biochemical_name',"class_metabolon","sub_class_metabolon","super_class_metabolon")
results_lvs<-merge(results_lvs, f_lvs_match, by = 'met', all = FALSE, sort = TRUE)

write.csv(results_lvs, "/udd/n2xwa/carb_met/results_lvs_20241205.csv")


#--------------------------------volcano plot ----------------------------------  

results_vol <- results_lvs
alpha_adjusted <- 0.05 / 293

# whole grain

results_vol <- transform(results_vol, color = ifelse(wgrain_intake == "positive", "#E41A1C",
                                         ifelse(wgrain_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(wgrain_intake == "positive" | wgrain_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_wgrain_adj < alpha_adjusted )

w <- ggplot(results_vol, 
            aes(x = beta_wgrain_adj, y = -log10(p_wgrain_adj), 
                color = color, label = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha_adjusted), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") + 
  geom_text_repel(data = significant_metabolites, aes(x = beta_wgrain_adj, y = -log10(p_wgrain_adj), label = metabolite_name), 
                  color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 3, 
                  hjust = 0.5, 
                  force = 3,
                  force_pull = 2,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf,
                  fontface = "bold") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold")) +
  labs(title = "Whole Grain Consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "Negative"))

ggsave("/udd/n2xwa/carb_met/volcano_wgrain.png", plot = w, width = 8, height = 8, dpi = 1200)


# select top ten: total carb
results_vol <- transform(results_vol, color = ifelse(carbo_intake == "positive", "#E41A1C",
                                                     ifelse(carbo_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(carbo_intake == "positive" | carbo_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_carb_adj < alpha_adjusted )

p2 = ggplot(results_vol, aes(x = beta_carb_adj, y = -log10(p_carb_adj), 
                                  color = color, label = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha_adjusted), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") + 
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  labs(title = "Total carbohydrate consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "negative"))  

# select top10
Up <- filter(results_vol, carbo_intake == 'positive') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_carb_adj))
Down <- filter(results_vol, carbo_intake == 'negative') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_carb_adj))

nudge_x_up = 0.01 - Up$beta_carb_adj 
nudge_x_down = -0.01 - Down$beta_carb_adj

p3 <- p2 + 
  geom_point(data = Up,aes(x = beta_carb_adj, y = -log10(p_carb_adj)),
             color = '#EB4232', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = beta_carb_adj, y = -log10(p_carb_adj), label = metabolite_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 3, 
                  nudge_x = nudge_x_up, 
                  direction = "y", 
                  hjust = 0, 
                  force = 3,
                  force_pull = 2,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf,
                  fontface = "bold")
p4 <- p3 + 
  geom_point(data = Down,aes(x = beta_carb_adj, y = -log10(p_carb_adj)),
             color = '#2DB2EB', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Down,aes(x = beta_carb_adj, y = -log10(p_carb_adj), label = metabolite_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 3, 
                  nudge_x = nudge_x_down, 
                  direction = "y", 
                  hjust = 0.5, 
                  force = 3,
                  force_pull = 2,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf,
                  fontface = "bold")

p4
ggsave("/udd/n2xwa/carb_met/volcano_carb.png", plot = p4, width = 8, height = 8, dpi = 1200)


# select top ten: added sugar

# select top ten
results_vol <- transform(results_vol, color = ifelse(addsugar_intake == "positive", "#E41A1C",
                                                     ifelse(addsugar_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(addsugar_intake == "positive" | addsugar_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_addsugar_adj < alpha_adjusted )

p2 = ggplot(results_vol, aes(x = beta_addsugar_adj, y = -log10(p_addsugar_adj), 
                             color = color, label = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha_adjusted), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") + 
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  labs(title = "Added sugar consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "negative"))  

# select top10
Up <- filter(results_vol, addsugar_intake == 'positive') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_addsugar_adj))
Down <- filter(results_vol, addsugar_intake == 'negative') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_addsugar_adj))

nudge_x_up = 0.005 - Up$beta_addsugar_adj 
nudge_x_down = -0.005 - Down$beta_addsugar_adj

p3 <- p2 + 
  geom_point(data = Up,aes(x = beta_addsugar_adj, y = -log10(p_addsugar_adj)),
             color = '#EB4232', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = beta_addsugar_adj, y = -log10(p_addsugar_adj), label = metabolite_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 3, 
                  nudge_x = nudge_x_up, 
                  direction = "y", 
                  hjust = 0, 
                  force = 3,
                  force_pull = 2,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf,
                  fontface = "bold")
p4 <- p3 + 
  geom_point(data = Down,aes(x = beta_addsugar_adj, y = -log10(p_addsugar_adj)),
             color = '#2DB2EB', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Down,aes(x = beta_addsugar_adj, y = -log10(p_addsugar_adj), label = metabolite_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 3, 
                  nudge_x = nudge_x_down, 
                  direction = "y", 
                  hjust = 0.5, 
                  force = 3,
                  force_pull = 2,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf,
                  fontface = "bold")

p4

ggsave("/udd/n2xwa/carb_met/volcano_addsugar.png", plot = p4, width = 8, height = 8, dpi = 1200)




#--------------------------------ploar plot ----------------------------------  

library(dplyr)
library(data.table)
library(tidyverse)
df <- results_vol
df
# df <- selected_metabolites
setDT(df)

# df <- df[, list(ID = metabolite_name, class =   class_metabolon, beta = beta_carb_adj, P =p_carb_adj )]

df <- subset(df, p_carb_adj < 0.05/293)
df <- df[, list(individual = metabolite_name, group =   class_metabolon, value = -log10(p_carb_adj)*sign(beta_carb_adj) )]
df
df$group <- ifelse(is.na(df$group), "NA", df$group)
# df <- df[order(group)]
# df$group <- factor(df$group, levels = unique(df$group))

df$group_label <- factor(df$group, levels = unique(df$group),labels = LETTERS[1:length(unique(df$group))])

df$group <- paste0(df$group_label,": ",  df$group)
df$group <- factor(df$group, levels = unique(df$group))

LETTERS[1:length(unique(df$group))]

library(tidyverse)
# Create dataset
data <- df

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

label_data

# prepare a data frame for base lines
base_data <-  data %>% 
  group_by(group_label) %>% dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

base_data


# prepare a data frame for grid (scales)
grid_data <- base_data

grid_data$end <- grid_data$end[ c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1

grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
grid_data


# Make the plot
alpha_c <-  0.2
label_data$individual <- sapply(label_data$individual, function(x) {
  # Use strsplit() and unlist() to split by "or" and take the first part
  sub("or.*", "", x)
})
# color = factor(value > 0)
polar <- ggplot(data, aes(x = as.factor(id), y = value, fill = group)) +
  geom_bar(aes(x = as.factor(id), y = value, fill = group), stat = "identity", alpha = 0.5) +
  
  # Add lines for 100/75/50/25 values
  geom_segment(data = grid_data, aes(x = end, y = 80, xend = start, yend = 80), 
               colour = "grey", alpha = alpha_c, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data, aes(x = end, y = 60, xend = start, yend = 60), 
               colour = "grey", alpha = alpha_c, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data, aes(x = end, y = 40, xend = start, yend = 40), 
               colour = "grey", alpha = alpha_c, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data, aes(x = end, y = 20, xend = start, yend = 20), 
               colour = "grey", alpha = alpha_c, size = 0.3, inherit.aes = FALSE) +
  
  # Annotate value lines
  annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), 
           label = c("20", "40", "60", "80"), color = "grey", size = 3, 
           angle = 0, fontface = "bold", hjust = 1) +
  scale_color_discrete(guide = "none") +
  
  ylim(-100, 120) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm"),
    text = element_text(face = "bold") 
  ) +
  coord_polar() +
  
  # Add labels
  geom_text(data = label_data, aes(x = id, y = value + 25, label = individual, hjust = hjust), 
            color = "black", fontface = "bold", alpha = 0.8, size = 3.7,
            angle = label_data$angle, inherit.aes = FALSE) +
  
  # Add base line information
  geom_segment(data = base_data, aes(x = start, y = 0, xend = end, yend = 0), 
               colour = "black", alpha = 0.5, size = 0.2, inherit.aes = FALSE) +
  geom_text(data = base_data, aes(x = title, y = -18, label = group_label), 
            colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)

polar
ggsave("/udd/n2xwa/carb_met/polar_carb.png", plot = polar, width = 12, height = 8, dpi = 1200)

#-------------------------------------------------------------------------------

#                     Elastic net regression model (ENR)

# Codes from Fenglei Wang /udd/nhfwa/Completed/PDI_met/Diabetologia

#-------------------------------------------------------------------------------              


names(temp_lvs)
lvs_ela <- temp_lvs[, c(1:294,355,379:388)]
names(lvs_ela)

# Convert data table to data frame
lvs_ela <- as.data.frame(lvs_ela)

#-------------------------------------------------------------------------------

#                      Elastic net for total carbohydrate

#-------------------------------------------------------------------------------


### training set and test set (70/30)   
set.seed(1234)
n=floor(nrow(lvs_ela)*0.7)
train_ind=sample(seq_len(nrow(lvs_ela)), size = n)
train=lvs_ela[train_ind,]
test=lvs_ela[-train_ind,]

#-------------------------------------------------------------------------------

#                      PREDICTION Training (total carbo)
#                          LOO CV APPROACH

#-------------------------------------------------------------------------------

record_carbo = data.frame(accumulateid=NA,newid=NA,NoMetabs=NA,accumulate_cor=NA)
x=0

for (i in 1:dim(train)[1]) {
  
  Training = cv.glmnet(as.matrix(train[-i,c(2:294)]), train[-i,"per_carb"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training$lambda.1se
  
  cmin = coef(Training, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  
  #### apply
  
  if(i==1) {
    train[i,dim(train)[2]+1] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  } else if(i>1) {
    
    train[i,dim(train)[2]] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  }
  
  #### test statistics
  
  x = x+1
  
  record_carbo[x,"accumulateid"] = i
  record_carbo[x,"newid"] = train[i,"id"]
  record_carbo[x,"NoMetabs"] = dim(metabsmin)[1]-1
  
  if(i<3) {
    record_carbo[x,"accumulate_cor"] = NA
  } else if(i>=3) {
    
    tmp = train[1:i,]
    record_carbo[x,"accumulate_cor"] = cor(tmp[,"per_carb"],
                                           tmp[,dim(tmp)[2]], use="complete.obs")
  }
  
}

record_carbo


#-------------------------------------------------------------------------------

#                   PREDICTION TESTING (carbo)

#-------------------------------------------------------------------------------

record_carbo = data.frame(repNo=NA,NoMetabs=NA,carbo_cor=NA, mse=NA)
met_coef_carbo = data.frame(met=c("intercept",colnames(train[,c(2:294)])))

for (inrep in 1:100) { 
  # traning set on the testing sets
  #cross-validation
  Training_CV = cv.glmnet(as.matrix(train[,c(2:294)]), train[,"per_carb"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training_CV$lambda.1se
  
  Training_M = glmnet(as.matrix(train[,c(2:294)]), train[,"per_carb"], family="gaussian", alpha=0.5) #normal distribution
  cmin = coef(Training_M, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  met_coef_carbo[,dim(met_coef_carbo)[2]+1]=cmin[,1]
  names(met_coef_carbo)[dim(met_coef_carbo)[2]] = paste("inrep",inrep,sep='')
  
  #### apply
  #############
  
  test[,dim(test)[2]+1] = as.numeric(predict(Training_M, as.matrix(test[,c(2:294)]), type="response", s=lambda_min_10F)[,1])
  
  test[,dim(test)[2]] = test[,dim(test)[2]] - metabsmin$Test.min[1]
  names(test)[dim(test)[2]] = paste("carbo",inrep,sep='')
  
  
  #### test statistics
  ####################
  
  record_carbo[inrep,"repNo"] = paste("inrep",inrep,sep='')
  record_carbo[inrep,"NoMetabs"] = dim(metabsmin)[1]-1
  record_carbo[inrep,"carbo_cor"] = cor(test[,"per_carb"],test[,dim(test)[2]])
  record_carbo[inrep,"mse"] = min(Training_CV$cvm)
  
}

record_carbo
met_coef_carbo


### find the score for carbo
table(record_carbo$NoMetabs) #  highest 30 times for 36 metabolites
#   33 36 38 41 47 50 60 65 71   #
#    3 30 27 18  8  6  3  1  4   #


sorted_data <- record_carbo %>%
  filter(NoMetabs == 36) %>%
  arrange(mse)

# inrep83: 36 metabolites (highest times) with the smallest MSE

met_carbo=as.character(met_coef_carbo[-1, ] %>% filter(inrep83!=0) %>% pull(met)) 
met_carbo
carbo_coef=as.character(met_coef_carbo[-1, ] %>% filter(inrep83!=0) %>% pull(inrep83)) 
carbo_coef

# create signature file
carbo_score_testing=test[,c("id", "per_carb","carbo83")] #carbo83 should be consistent with (inrep83
colnames(carbo_score_testing)[3]="carbo_score"
carbo_score_training <- train[, c("id", "per_carb", grep("^V", names(train), value = TRUE))]
colnames(carbo_score_training)[3]="carbo_score"

identical(colnames(carbo_score_testing), colnames(carbo_score_training))
carbo_score_testing$set="Testing"
carbo_score_training$set="Training"

signature_carbo=rbind(carbo_score_testing, carbo_score_training) %>% arrange(id)
signature_carbo
write.csv(signature_carbo, "/udd/n2xwa/carb_met/signature_carbo.csv", na="")

# correlation between metabolite profile score and carbohydrate intake #
cor_pearson_train <- cor(carbo_score_training$carbo_score, carbo_score_training$per_carb, method = "pearson")
cor_spearman_train <- cor(carbo_score_training$carbo_score, carbo_score_training$per_carb, method = "spearman")
cor_pearson_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$per_carb, method = "pearson")
cor_spearman_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$per_carb, method = "spearman")

print(cor_pearson_train)
print(cor_spearman_train)
print(cor_pearson_test)
print(cor_spearman_test)


# save ENR results to "XXX.RData" and later load
carbo_ela <- list(test, train, record_carbo, met_coef_carbo, met_carbo, carbo_coef, Training_M)
save(carbo_ela, file = "/udd/n2xwa/carb_met/carbo_ela_293.RData")


