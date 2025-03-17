#-------------------------------------------------------------------------------
#  December 7, 2024
#  The purpose of this program is to examine the association between metabolites signature 
#  and T2D among NHS, NHSII, and HPFS 
#  Programmer: Xiaowen Wang /udd/n2xwa/carb_met/
#-------------------------------------------------------------------------------

library(survival)

#-------------------------------------------------------------------------------

# [1]. load blood cohort metabolites, outcomes, and covariates

#-------------------------------------------------------------------------------


load(file="/udd/n2xwa/metabolomics_first_13150_imputed.RData")
f.cohort <- feature_met1
p.cohort <- pheno_met1_imputed
a.cohort <- assay_met1_imputed

z.cohort=cbind(a.cohort,p.cohort) 
z.cohort=z.cohort %>%  mutate(id=str_sub(id, start = 1, end = 6),
                        newid=case_when(cohort=="nhs1" ~ as.numeric(id)+1000000,
                                        cohort=="nhs2" ~ as.numeric(id)+2000000,
                                        cohort=="hpfs" ~ as.numeric(id)+3000000))
names(z.cohort)


### OUtcome from cohort linked to those who have blood data #####
p_bld<-p.cohort %>% mutate(id=str_sub(id, start = 1, end = 6)) %>%  dplyr::select(id,blddate,cohort)
write.csv(p_bld, "/udd/n2xwa/lvs_met/dbld_all.csv", na="", row.names = F)
haven::write_sas(p_bld, "/udd/n2xwa/lvs_met/dbld_all.sas7bdat")
foreign::write.foreign(df = p_bld,
                       datafile = '/udd/n2xwa/lvs_met/dbld_all.txt',
                       codefile = '/udd/n2xwa/lvs_met/dbld_all.sas',
                       dataname = '/udd/n2xwa/lvs_met/', # Destination in SAS to save the data
                       package = 'SAS')

# OUtcome and cov from cohort
# n1_cov.csv,n2_cov.csv, h_cov.csv were created by n1_cov.sas, n2_cov.sas, h_cov.sas, respectively
n1_cov<-read.csv("/udd/n2xwa/lvs_met/cohort/n1_cov.csv", header = TRUE)
n2_cov<-read.csv("/udd/n2xwa/lvs_met/cohort/n2_cov.csv", header = TRUE)
h_cov<-read.csv("/udd/n2xwa/lvs_met/cohort/h_cov.csv", header = TRUE)

n1_cov<-n1_cov %>% mutate(newid=as.numeric(id)+1000000,bmi_young=as.numeric(bmi18))
n2_cov<-n2_cov %>% mutate(newid=as.numeric(id)+2000000,bmi_young=as.numeric(bmi18))
h_cov<-h_cov %>% mutate(newid=as.numeric(id)+3000000,bmi_young=as.numeric(bmi2186))


z_n1<-z.cohort %>% filter(cohort=='nhs1')
z_n2<-z.cohort %>% filter(cohort=='nhs2')
z_h <-z.cohort %>% filter(cohort=='hpfs')

z_n1<-merge(z_n1, n1_cov, by = 'newid',all = F, sort = TRUE) 
z_n2<-merge(z_n2, n2_cov, by = 'newid',all = F, sort = TRUE)
z_h<-merge(z_h, h_cov, by = 'newid',all = F, sort = TRUE)
z_h<-z_h %>% mutate(phms_bld=NA)


met=dput(names(z_n1[2:261]))

temp_n1<-z_n1 %>% dplyr::select("newid",all_of(met),"cohort","caco", "matchid", "fast","race","dbfh","age_bld", "bmi_bld", "smk_bld", "phms_bld", 
                                "alco_bld", "ahei_bld","aheinowgr_bld", "calor_bld", "act_bld", "hbp_bld", "chol_bld","bmi_young", 
                                "carbo_bld","prot_bld",  "pot_bld" , "leg_bld", "whg_bld", "ref_bld","sumref_bld","fruit_bld", "veg_bld", "sugar_bld",
                                 "mifh", "blddate","dtdxdb2","dtdxcvd","dtdxchd","dtdxstr","dtdth","cutoff","db2case","tdb2","cvdcase","tcvd",
                                "db_out", "cvd_out", "chd_out", "str_out", "dth_out", "dead_cvd", "dead_ca", "dead_db")
temp_n2<-z_n2 %>% dplyr::select("newid",all_of(met),"cohort","caco", "matchid", "fast","race","dbfh","age_bld", "bmi_bld", "smk_bld", "phms_bld", 
                                "alco_bld", "ahei_bld", "aheinowgr_bld", "calor_bld", "act_bld", "hbp_bld", "chol_bld","bmi_young", 
                                "carbo_bld","prot_bld",  "pot_bld" , "leg_bld", "whg_bld", "ref_bld","sumref_bld","fruit_bld", "veg_bld", "sugar_bld",
                                "mifh", "blddate","dtdxdb2","dtdxcvd","dtdxchd","dtdxstr","dtdth","cutoff", "db2case","tdb2","cvdcase","tcvd",
                                "db_out", "cvd_out", "chd_out", "str_out", "dth_out", "dead_cvd", "dead_ca", "dead_db")


temp_h<-z_h %>% dplyr::select("newid",all_of(met),"cohort","caco", "matchid", "fast","race","dbfh","age_bld", "bmi_bld", "smk_bld","phms_bld", 
                              "alco_bld", "ahei_bld","aheinowgr_bld", "calor_bld", "act_bld", "hbp_bld", "chol_bld","bmi_young", 
                              "carbo_bld","prot_bld",  "pot_bld" , "leg_bld", "whg_bld", "ref_bld","sumref_bld","fruit_bld", "veg_bld", "sugar_bld",
                              "mifh", "blddate","dtdxdb2","dtdxcvd","dtdxchd","dtdxstr","dtdth","cutoff", "db2case","tdb2","cvdcase","tcvd",
                              "db_out", "cvd_out", "chd_out", "str_out", "dth_out", "dead_cvd", "dead_ca", "dead_db")

temp_cohort=rbind(temp_n1,temp_n2,temp_h) # 11454 baseline data set


#-------------------------------------------------------------------------------
#   [2].  carb predicted score and T2D:
#-------------------------------------------------------------------------------
  # (1)  Load ENR results, respectively:

# 1. carbo_ela_293.RData (Total)
# 2. addsugar_ela_293.RData
# 3. wgrain_ela_293.RData
# 4. rgrain_ela_293.RData
# 5. legume_ela_293.RData
# 6. potato_ela_293.RData
# 7. fruit_ela_293.RData
# 8. vegetable_ela_293.RData
 
#-------------------------------------------------------------------------------
load(file = "/udd/n2xwa/carb_met/carbo_ela_293.RData")
list(carbo_ela)
record_carbo <- carbo_ela[[3]]
Training_M <- carbo_ela[[7]]
met_carbo <- carbo_ela[[5]]
carbo_coef <- carbo_ela[[6]]
met_coef_carbo <-carbo_ela[[4]]

#-------------------------------------------------------------------------------

# (2)  using parametric method to calculate predicted score 

#-------------------------------------------------------------------------------
# check the metabolites in X and B
df_beta <- data.frame(metab_name =met_carbo,coef= as.numeric(carbo_coef))
df_beta

df_data <- temp_cohort

gplots::venn(list(names(df_data),
                  df_beta$metab_name
))

#  metabolites overlapping

v_overlap <- intersect(names(df_data ),df_beta$metab_name)

df_beta <- df_beta[df_beta$metab_name %in% v_overlap, ]
df_data
df_data <- df_data[,c("newid", v_overlap) ]


table(names(df_data)[-1] == df_beta$metab_name)


# matrix multiply

df_data_Y <- as.matrix(df_data[, -1]) %*% as.matrix(df_beta$coef)
df_data_Y

df_data_Y <- data.frame(newid = df_data$newid, y_predict = df_data_Y)
temp_cohort<-merge(temp_cohort, df_data_Y, by = 'newid', all = FALSE, sort = TRUE)
names(temp_cohort)
temp_cohort$y_predictSD <- scale(temp_cohort$y_predict)


#-------------------------------------------------------------------------------

# (3)  predicted carb met score and T2D: cox models

#-------------------------------------------------------------------------------


str(temp_cohort$calor_bld)
temp_cohort$calor_bld <- as.numeric(temp_cohort$calor_bld)
temp_cohort<-temp_cohort %>% group_by(cohort) %>% mutate(white=case_when(race=='White' ~ 1,T ~ 0),
                                                         fasting=case_when(fast=='nonfasting' ~ 0,T ~ 1),
                                                         smkstatus=case_when(smk_bld==3 ~ 'current',
                                                                             smk_bld==2 ~ 'past',
                                                                             T ~ 'never'),
                                                         neversmoking=case_when(smk_bld>1 ~ 'no',
                                                                                T ~ 'yes'),
                                                         sex=case_when(cohort=='hpfs' ~ 'male',
                                                                       T ~ 'female'),
                                                         bmi_group=case_when(bmi_bld<25 ~ 'normal',
                                                                             bmi_bld<30 ~ 'overweight',
                                                                             bmi_bld>=30 ~ 'obese',
                                                                             T ~ 'NA'),
                                                         prot_bld = tidyr::replace_na(prot_bld, median(prot_bld, na.rm = TRUE)),
                                                         carbo_bld = tidyr::replace_na(carbo_bld, median(carbo_bld, na.rm = TRUE)),
                                                         pot_bld = tidyr::replace_na(pot_bld, median(pot_bld, na.rm = TRUE)),
                                                         whg_bld = tidyr::replace_na(whg_bld, median(whg_bld, na.rm = TRUE)),
                                                         leg_bld = tidyr::replace_na(leg_bld, median(leg_bld, na.rm = TRUE)),
                                                         ref_bld = tidyr::replace_na(ref_bld, median(ref_bld, na.rm = TRUE)),
                                                         sugar_bld = tidyr::replace_na(sugar_bld, median(sugar_bld, na.rm = TRUE)),
                                                         sumref_bld = tidyr::replace_na(sumref_bld, median(sumref_bld, na.rm = TRUE)),
                                                         fruit_bld = tidyr::replace_na(fruit_bld, median(fruit_bld, na.rm = TRUE)),
                                                         veg_bld = tidyr::replace_na(veg_bld, median(veg_bld, na.rm = TRUE)),
                                                         act_bld = tidyr::replace_na(act_bld, median(act_bld, na.rm = TRUE)),
                                                         alco_bld = tidyr::replace_na(alco_bld, median(alco_bld, na.rm = TRUE)),
                                                         ahei_bld = tidyr::replace_na(ahei_bld, median(ahei_bld, na.rm = TRUE)),
                                                         aheinowgr_bld = tidyr::replace_na(aheinowgr_bld, median(aheinowgr_bld, na.rm = TRUE)),
                                                         calor_bld = tidyr::replace_na(calor_bld, median(calor_bld, na.rm = TRUE)),
                                                         alco_g=statar::xtile(alco_bld, n = 5),
                                                         act_g=statar::xtile(act_bld, n = 5),
                                                         ahei_g=statar::xtile(ahei_bld, n = 5),
                                                         aheinowgr_g=statar::xtile(aheinowgr_bld, n = 5),
                                                         prot_g=statar::xtile(prot_bld, n = 5),
                                                         calor_g=statar::xtile(calor_bld, n = 5)) %>% as.data.frame()

temp_cohort<-temp_cohort %>% mutate(dtdxdb2=as.numeric(dtdxdb2),
                                    dtdxcvd=as.numeric(dtdxcvd),
                                    dtdxchd=as.numeric(dtdxchd),
                                    dtdxstr=as.numeric(dtdxstr),
                                    dtdth=as.numeric(dtdth),
                                    cutoff=as.numeric(cutoff),
                                    blddate=as.numeric(blddate),
                                    dtdxdb2=case_when(is.na(dtdxdb2) ~ 9999, T ~ dtdxdb2),
                                    dtdth=case_when(is.na(dtdth) ~ 9999, T ~ dtdth),
                                    dtdxcvd=case_when(is.na(dtdxcvd) ~ 9999, T ~ dtdxcvd))  


# temp_cohort<-temp_cohort %>% mutate(db_time=case_when(db_out==1 ~ dtdxdb2-blddate,T ~ cutoff-blddate),
#                                    dth_time=case_when(dth_out==1 ~ dtdth-blddate,T ~ cutoff-blddate),
#                                    cvd_time=case_when(cvd_out==1 ~ dtdxcvd-blddate,T ~ cutoff-blddate),
#                                    chd_time=case_when(chd_out==1 ~ dtdxchd-blddate,T ~ cutoff-blddate),
#                                    str_time=case_when(str_out==1 ~ dtdxstr-blddate,T ~ cutoff-blddate))

temp_cohort<-temp_cohort %>% filter(!duplicated(newid)) # id no duplicate

table(temp_cohort$db2case,temp_cohort$cohort)
table(temp_cohort$cvdcase,temp_cohort$cohort)
# temp_cohort <- temp_cohort %>% mutate(across(c(carbo_bld, prot_bld, pot_bld, leg_bld, whg_bld, ref_bld),
#                                             ~ . * 4 / calor_bld, .names = "{.col}_energy"))


## Metabolic score and carb in cohort <<<<< box-plot
library(Hmisc)
temp_cohort$carbo_bld_5 <- cut(temp_cohort$carbo_bld, breaks = 5)
Rho <- spearman(temp_cohort$carbo_bld_5, temp_cohort$y_predictSD)
box <- ggplot(temp_cohort, aes(x = carbo_bld_5, y = y_predictSD)) +
  geom_boxplot(fill = "#add8e6") +
  xlab("Quintiles of total carbohydrate intake") +
  ylab("Metabolic score") +
  # 添加Spearman相关系数和P值（文本加粗）
  geom_text(x = 5, y = max(temp_cohort$y_predictSD), 
            label = paste("Spearman Rho =", round(Rho, 3), "\nP < 0.0001"),
            hjust = 1, vjust = 1, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold"), # 坐标轴标题加粗
    axis.text = element_text(face = "bold"),  # 坐标轴文字加粗
    #plot.title = element_text(face = "bold"),  # 图标题加粗 (如果需要)
  )

ggsave("/udd/n2xwa/carb_met/box_carb.png", plot = box, width = 8, height = 8, dpi = 1200)

### Cox regression
fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld, data=temp_cohort)
fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld  +
            dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
            factor(act_g) + factor(alco_g)+ factor(hbp_bld) + 
            factor(chol_bld) , data=temp_cohort)

fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld  + bmi_young +
            dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
            factor(act_g) + factor(alco_g)+ factor(hbp_bld) + factor(prot_g)+
            factor(chol_bld) , data=temp_cohort)


R<-data.frame(variable='model 1',HR=summary(fit)$conf.int[,1],LL=summary(fit)$conf.int[,3],
              UL=summary(fit)$conf.int[,4],P=summary(fit)$coefficients[,5])


R<-R %>% mutate(HRCI=paste0(round(HR,3),"(",round(LL,3),"-",round(UL,3),")"),Pvalue=round(P,3))


R
# total carb
#                   variable        HR        LL        UL            P               HRCI Pvalue
#y_predictSD         model 1 1.0674137 1.0145871 1.1229909 1.176313e-02 1.067(1.015-1.123)  0.012

### quintile HR 
#  y_predictSD quintile (1, 2, 3, 4, 5)
temp_cohort$y_predictSD_quintile <- ntile(temp_cohort$y_predictSD, 5)
temp_cohort$y_predictSD_quintile <- factor(temp_cohort$y_predictSD_quintile, levels = 1:5)
fit1_q <- coxph(Surv(tdb2, db2case) ~ y_predictSD_quintile + cohort + age_bld, data = temp_cohort)
fit2_q <- coxph(Surv(tdb2, db2case) ~ y_predictSD_quintile + cohort + age_bld + 
                  dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                  factor(act_g) + factor(alco_g) + factor(hbp_bld) + factor(chol_bld), 
                data = temp_cohort)
dim(summary(fit1_q)$conf.int)
dim(summary(fit1_q)$coefficients)
rownames(summary(fit1_q)$coefficients)

quintile_rows <- grep("y_predictSD_quintile", rownames(summary(fit1_q)$conf.int))
fit1_confint <- summary(fit1_q)$conf.int[quintile_rows, ]
fit1_coeff <- summary(fit1_q)$coefficients[quintile_rows, ]

R_fit1 <- data.frame(
  variable = paste0("Q", 2:5),  #  Q2-Q5
  HR = fit1_confint[, 1],       #  HR
  LL = fit1_confint[, 3],       # 
  UL = fit1_confint[, 4],       # 
  P = fit1_coeff[, 5]           # P 
)

R_fit1 <- R_fit1 %>%
  mutate(HRCI = paste0(round(HR, 2), "(", round(LL, 2), "-", round(UL, 2), ")"),
         Pvalue = round(P, 2))


quintile_rows_fit2 <- grep("y_predictSD_quintile", rownames(summary(fit2_q)$conf.int))
fit2_confint <- summary(fit2_q)$conf.int[quintile_rows_fit2, ]
fit2_coeff <- summary(fit2_q)$coefficients[quintile_rows_fit2, ]
R_fit2 <- data.frame(
  variable = paste0("Q", 2:5),  #  Q2-Q5 
  HR = fit2_confint[, 1],       #  HR
  LL = fit2_confint[, 3],        
  UL = fit2_confint[, 4],       
  P = fit2_coeff[, 5]           #  P 
)

R_fit2 <- R_fit2 %>%
  mutate(HRCI = paste0(round(HR, 2), "(", round(LL, 2), "-", round(UL, 2), ")"),
         Pvalue = round(P, 2))

R_combined <- rbind(R_fit1, R_fit2)
R_combined

# variable       HR        LL       UL            P            HRCI Pvalue
# y_predictSD_quintile2        Q2 1.208219 1.0327374 1.413519 1.816324e-02 1.21(1.03-1.41)   0.02
# y_predictSD_quintile3        Q3 1.239500 1.0604853 1.448734 6.978143e-03 1.24(1.06-1.45)   0.01
# y_predictSD_quintile4        Q4 1.211324 1.0360616 1.416235 1.620617e-02 1.21(1.04-1.42)   0.02
# y_predictSD_quintile5        Q5 1.389647 1.1930068 1.618699 2.368302e-05 1.39(1.19-1.62)   0.00
# y_predictSD_quintile21       Q2 1.105876 0.9442055 1.295229 2.120284e-01  1.11(0.94-1.3)   0.21
# y_predictSD_quintile31       Q3 1.111491 0.9496752 1.300879 1.879255e-01  1.11(0.95-1.3)   0.19
# y_predictSD_quintile41       Q4 1.050663 0.8968206 1.230896 5.406546e-01  1.05(0.9-1.23)   0.54
# y_predictSD_quintile51       Q5 1.168993 0.9995577 1.367149 5.065085e-02    1.17(1-1.37)   0.05


# P for trend
quintile_medians <- temp_cohort %>%
  group_by(y_predictSD_quintile) %>%
  summarise(median_value = median(y_predictSD, na.rm = TRUE))
temp_cohort <- temp_cohort %>%
  left_join(quintile_medians, by = c("y_predictSD_quintile" = "y_predictSD_quintile"))


fit1_trend <- coxph(Surv(tdb2, db2case) ~ median_value + cohort + age_bld, data = temp_cohort)
fit2_trend <- coxph(Surv(tdb2, db2case) ~ median_value + cohort + age_bld + 
                      dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                      factor(act_g) + factor(alco_g) + factor(hbp_bld) + factor(chol_bld), 
                    data = temp_cohort)
p_for_trend_fit1 <- summary(fit1_trend)$coefficients["median_value", "Pr(>|z|)"]
p_for_trend_fit2 <- summary(fit2_trend)$coefficients["median_value", "Pr(>|z|)"]
p_for_trend_fit1
p_for_trend_fit2




#-------------------------------------------------------------------------------
#   [3].  Repeat the above process: carb predicted score and T2D:
#-------------------------------------------------------------------------------
# (1)  Load ENR results:

# 1. wgrain_ela_293.RData

#-------------------------------------------------------------------------------
load(file = "/udd/n2xwa/carb_met/cohort/wgrain_ela_293.RData")
list(carbo_ela)
record_carbo <- carbo_ela[[3]]
Training_M <- carbo_ela[[7]]
met_carbo <- carbo_ela[[5]]
carbo_coef <- carbo_ela[[6]]
met_coef_carbo <-carbo_ela[[4]]

temp_cohort=rbind(temp_n1,temp_n2,temp_h) # 11454 baseline data set

#-------------------------------------------------------------------------------

# (2)  using parametric method to calculate predicted score 

#-------------------------------------------------------------------------------
# check the metabolites in X and B
df_beta <- data.frame(metab_name =met_carbo,coef= as.numeric(carbo_coef))
df_beta

df_data <- temp_cohort

gplots::venn(list(names(df_data),
                  df_beta$metab_name
))

#  metabolites overlapping

v_overlap <- intersect(names(df_data ),df_beta$metab_name)

df_beta <- df_beta[df_beta$metab_name %in% v_overlap, ]
df_data
df_data <- df_data[,c("newid", v_overlap) ]


table(names(df_data)[-1] == df_beta$metab_name)


# matrix multiply

df_data_Y <- as.matrix(df_data[, -1]) %*% as.matrix(df_beta$coef)
df_data_Y

df_data_Y <- data.frame(newid = df_data$newid, y_predict = df_data_Y)
temp_cohort<-merge(temp_cohort, df_data_Y, by = 'newid', all = FALSE, sort = TRUE)
names(temp_cohort)
temp_cohort$y_predictSD <- scale(temp_cohort$y_predict)


#-------------------------------------------------------------------------------

# (3)  predicted carb met score and T2D: cox models

#-------------------------------------------------------------------------------


str(temp_cohort$calor_bld)
temp_cohort$calor_bld <- as.numeric(temp_cohort$calor_bld)
temp_cohort<-temp_cohort %>% group_by(cohort) %>% mutate(white=case_when(race=='White' ~ 1,T ~ 0),
                                                         fasting=case_when(fast=='nonfasting' ~ 0,T ~ 1),
                                                         smkstatus=case_when(smk_bld==3 ~ 'current',
                                                                             smk_bld==2 ~ 'past',
                                                                             T ~ 'never'),
                                                         neversmoking=case_when(smk_bld>1 ~ 'no',
                                                                                T ~ 'yes'),
                                                         sex=case_when(cohort=='hpfs' ~ 'male',
                                                                       T ~ 'female'),
                                                         bmi_group=case_when(bmi_bld<25 ~ 'normal',
                                                                             bmi_bld<30 ~ 'overweight',
                                                                             bmi_bld>=30 ~ 'obese',
                                                                             T ~ 'NA'),
                                                         prot_bld = tidyr::replace_na(prot_bld, median(prot_bld, na.rm = TRUE)),
                                                         carbo_bld = tidyr::replace_na(carbo_bld, median(carbo_bld, na.rm = TRUE)),
                                                         pot_bld = tidyr::replace_na(pot_bld, median(pot_bld, na.rm = TRUE)),
                                                         whg_bld = tidyr::replace_na(whg_bld, median(whg_bld, na.rm = TRUE)),
                                                         leg_bld = tidyr::replace_na(leg_bld, median(leg_bld, na.rm = TRUE)),
                                                         ref_bld = tidyr::replace_na(ref_bld, median(ref_bld, na.rm = TRUE)),
                                                         sugar_bld = tidyr::replace_na(sugar_bld, median(sugar_bld, na.rm = TRUE)),
                                                         sumref_bld = tidyr::replace_na(sumref_bld, median(sumref_bld, na.rm = TRUE)),
                                                         fruit_bld = tidyr::replace_na(fruit_bld, median(fruit_bld, na.rm = TRUE)),
                                                         veg_bld = tidyr::replace_na(veg_bld, median(veg_bld, na.rm = TRUE)),
                                                         act_bld = tidyr::replace_na(act_bld, median(act_bld, na.rm = TRUE)),
                                                         alco_bld = tidyr::replace_na(alco_bld, median(alco_bld, na.rm = TRUE)),
                                                         ahei_bld = tidyr::replace_na(ahei_bld, median(ahei_bld, na.rm = TRUE)),
                                                         aheinowgr_bld = tidyr::replace_na(aheinowgr_bld, median(aheinowgr_bld, na.rm = TRUE)),
                                                         calor_bld = tidyr::replace_na(calor_bld, median(calor_bld, na.rm = TRUE)),
                                                         alco_g=statar::xtile(alco_bld, n = 5),
                                                         act_g=statar::xtile(act_bld, n = 5),
                                                         ahei_g=statar::xtile(ahei_bld, n = 5),
                                                         aheinowgr_g=statar::xtile(aheinowgr_bld, n = 5),
                                                         prot_g=statar::xtile(prot_bld, n = 5),
                                                         calor_g=statar::xtile(calor_bld, n = 5)) %>% as.data.frame()

temp_cohort<-temp_cohort %>% mutate(dtdxdb2=as.numeric(dtdxdb2),
                                    dtdxcvd=as.numeric(dtdxcvd),
                                    dtdxchd=as.numeric(dtdxchd),
                                    dtdxstr=as.numeric(dtdxstr),
                                    dtdth=as.numeric(dtdth),
                                    cutoff=as.numeric(cutoff),
                                    blddate=as.numeric(blddate),
                                    dtdxdb2=case_when(is.na(dtdxdb2) ~ 9999, T ~ dtdxdb2),
                                    dtdth=case_when(is.na(dtdth) ~ 9999, T ~ dtdth),
                                    dtdxcvd=case_when(is.na(dtdxcvd) ~ 9999, T ~ dtdxcvd))  


# temp_cohort<-temp_cohort %>% mutate(db_time=case_when(db_out==1 ~ dtdxdb2-blddate,T ~ cutoff-blddate),
#                                    dth_time=case_when(dth_out==1 ~ dtdth-blddate,T ~ cutoff-blddate),
#                                    cvd_time=case_when(cvd_out==1 ~ dtdxcvd-blddate,T ~ cutoff-blddate),
#                                    chd_time=case_when(chd_out==1 ~ dtdxchd-blddate,T ~ cutoff-blddate),
#                                    str_time=case_when(str_out==1 ~ dtdxstr-blddate,T ~ cutoff-blddate))

temp_cohort<-temp_cohort %>% filter(!duplicated(newid)) # id no duplicate

table(temp_cohort$db2case,temp_cohort$cohort)
table(temp_cohort$cvdcase,temp_cohort$cohort)
# temp_cohort <- temp_cohort %>% mutate(across(c(carbo_bld, prot_bld, pot_bld, leg_bld, whg_bld, ref_bld),
#                                             ~ . * 4 / calor_bld, .names = "{.col}_energy"))


### Cox regression
fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld, data=temp_cohort)
fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld  +
            dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
            factor(act_g) + factor(alco_g)+ factor(hbp_bld) + 
            factor(chol_bld) , data=temp_cohort)

fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld  + bmi_young +
            dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
            factor(act_g) + factor(alco_g)+ factor(hbp_bld) + 
            factor(chol_bld) , data=temp_cohort)


R<-data.frame(variable='model 1',HR=summary(fit)$conf.int[,1],LL=summary(fit)$conf.int[,3],
              UL=summary(fit)$conf.int[,4],P=summary(fit)$coefficients[,5])


R<-R %>% mutate(HRCI=paste0(round(HR,3),"(",round(LL,3),"-",round(UL,3),")"),Pvalue=round(P,3))


R
# total carb
#                   variable        HR        LL        UL            P               HRCI Pvalue
#y_predictSD        model 1 0.7293798 0.6941788 0.7663659 7.140960e-36 0.729(0.694-0.766)  0.000

### quintile HR 
#  y_predictSD quintile (1, 2, 3, 4, 5)
temp_cohort$y_predictSD_quintile <- ntile(temp_cohort$y_predictSD, 5)
temp_cohort$y_predictSD_quintile <- factor(temp_cohort$y_predictSD_quintile, levels = 1:5)
fit1_q <- coxph(Surv(tdb2, db2case) ~ y_predictSD_quintile + cohort + age_bld, data = temp_cohort)
fit2_q <- coxph(Surv(tdb2, db2case) ~ y_predictSD_quintile + cohort + age_bld + 
                  dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                  factor(act_g) + factor(alco_g) + factor(hbp_bld) + factor(chol_bld), 
                data = temp_cohort)
dim(summary(fit1_q)$conf.int)
dim(summary(fit1_q)$coefficients)
rownames(summary(fit1_q)$coefficients)

quintile_rows <- grep("y_predictSD_quintile", rownames(summary(fit1_q)$conf.int))


fit1_confint <- summary(fit1_q)$conf.int[quintile_rows, ]
fit1_coeff <- summary(fit1_q)$coefficients[quintile_rows, ]

R_fit1 <- data.frame(
  variable = paste0("Q", 2:5),  #  Q2-Q5
  HR = fit1_confint[, 1],       #  HR
  LL = fit1_confint[, 3],       # 
  UL = fit1_confint[, 4],       # 
  P = fit1_coeff[, 5]           # P 
)

R_fit1 <- R_fit1 %>%
  mutate(HRCI = paste0(round(HR, 2), "(", round(LL, 2), "-", round(UL, 2), ")"),
         Pvalue = round(P, 3))


quintile_rows_fit2 <- grep("y_predictSD_quintile", rownames(summary(fit2_q)$conf.int))
fit2_confint <- summary(fit2_q)$conf.int[quintile_rows_fit2, ]
fit2_coeff <- summary(fit2_q)$coefficients[quintile_rows_fit2, ]
R_fit2 <- data.frame(
  variable = paste0("Q", 2:5),  #  Q2-Q5 
  HR = fit2_confint[, 1],       #  HR
  LL = fit2_confint[, 3],        
  UL = fit2_confint[, 4],       
  P = fit2_coeff[, 5]           #  P 
)

R_fit2 <- R_fit2 %>%
  mutate(HRCI = paste0(round(HR, 2), "(", round(LL, 2), "-", round(UL, 2), ")"),
         Pvalue = round(P, 3))

R_combined <- rbind(R_fit1, R_fit2)
R_combined

# variable       HR        LL       UL            P            HRCI Pvalue
# y_predictSD_quintile2        Q2 0.7209743 0.6319380 0.8225553 1.147237e-06 0.72(0.63-0.82)      0
# y_predictSD_quintile3        Q3 0.5581952 0.4845733 0.6430027 6.508966e-16 0.56(0.48-0.64)      0
# y_predictSD_quintile4        Q4 0.4993775 0.4316820 0.5776888 9.387680e-21  0.5(0.43-0.58)      0
# y_predictSD_quintile5        Q5 0.3818919 0.3257136 0.4477598 1.976610e-32 0.38(0.33-0.45)      0
# y_predictSD_quintile21       Q2 0.7299973 0.6394990 0.8333025 3.156217e-06 0.73(0.64-0.83)      0
# y_predictSD_quintile31       Q3 0.6016730 0.5218001 0.6937723 2.727035e-12  0.6(0.52-0.69)      0
# y_predictSD_quintile41       Q4 0.5411027 0.4671557 0.6267550 2.581126e-16 0.54(0.47-0.63)      0
# y_predictSD_quintile51       Q5 0.4202937 0.3576930 0.4938504 6.026204e-26 0.42(0.36-0.49)      0


# P for trend
quintile_medians <- temp_cohort %>%
  group_by(y_predictSD_quintile) %>%
  summarise(median_value = median(y_predictSD, na.rm = TRUE))
temp_cohort <- temp_cohort %>%
  left_join(quintile_medians, by = c("y_predictSD_quintile" = "y_predictSD_quintile"))


fit1_trend <- coxph(Surv(tdb2, db2case) ~ median_value + cohort + age_bld, data = temp_cohort)
fit2_trend <- coxph(Surv(tdb2, db2case) ~ median_value + cohort + age_bld + 
                      dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                      factor(act_g) + factor(alco_g) + factor(hbp_bld) + factor(chol_bld), 
                    data = temp_cohort)
p_for_trend_fit1 <- summary(fit1_trend)$coefficients["median_value", "Pr(>|z|)"]
p_for_trend_fit2 <- summary(fit2_trend)$coefficients["median_value", "Pr(>|z|)"]
p_for_trend_fit1
p_for_trend_fit2


#-------------------------------------------------------------------------------
#   [4].  Repeat the above process: carb predicted score and T2D:
#-------------------------------------------------------------------------------
# (1)  Load ENR results:

# 1. addsugar_ela_293.RData

#-------------------------------------------------------------------------------
load(file = "/udd/n2xwa/carb_met/cohort/addsugar_ela_293.RData")
list(carbo_ela)
record_carbo <- carbo_ela[[3]]
Training_M <- carbo_ela[[7]]
met_carbo <- carbo_ela[[5]]
carbo_coef <- carbo_ela[[6]]
met_coef_carbo <-carbo_ela[[4]]

temp_cohort=rbind(temp_n1,temp_n2,temp_h) # 11454 baseline data set

#-------------------------------------------------------------------------------

# (2)  using parametric method to calculate predicted score 

#-------------------------------------------------------------------------------
# check the metabolites in X and B

df_beta <- data.frame(metab_name =met_carbo,coef= as.numeric(carbo_coef))
df_beta

df_data <- temp_cohort

gplots::venn(list(names(df_data),
                  df_beta$metab_name
))

#  metabolites overlapping

v_overlap <- intersect(names(df_data ),df_beta$metab_name)

df_beta <- df_beta[df_beta$metab_name %in% v_overlap, ]
df_data
df_data <- df_data[,c("newid", v_overlap) ]


table(names(df_data)[-1] == df_beta$metab_name)


# matrix multiply

df_data_Y <- as.matrix(df_data[, -1]) %*% as.matrix(df_beta$coef)
df_data_Y

df_data_Y <- data.frame(newid = df_data$newid, y_predict = df_data_Y)
temp_cohort<-merge(temp_cohort, df_data_Y, by = 'newid', all = FALSE, sort = TRUE)
names(temp_cohort)
temp_cohort$y_predictSD <- scale(temp_cohort$y_predict)


#-------------------------------------------------------------------------------

# (3)  predicted carb met score and T2D: cox models

#-------------------------------------------------------------------------------


str(temp_cohort$calor_bld)
temp_cohort$calor_bld <- as.numeric(temp_cohort$calor_bld)
temp_cohort<-temp_cohort %>% group_by(cohort) %>% mutate(white=case_when(race=='White' ~ 1,T ~ 0),
                                                         fasting=case_when(fast=='nonfasting' ~ 0,T ~ 1),
                                                         smkstatus=case_when(smk_bld==3 ~ 'current',
                                                                             smk_bld==2 ~ 'past',
                                                                             T ~ 'never'),
                                                         neversmoking=case_when(smk_bld>1 ~ 'no',
                                                                                T ~ 'yes'),
                                                         sex=case_when(cohort=='hpfs' ~ 'male',
                                                                       T ~ 'female'),
                                                         bmi_group=case_when(bmi_bld<25 ~ 'normal',
                                                                             bmi_bld<30 ~ 'overweight',
                                                                             bmi_bld>=30 ~ 'obese',
                                                                             T ~ 'NA'),
                                                         prot_bld = tidyr::replace_na(prot_bld, median(prot_bld, na.rm = TRUE)),
                                                         carbo_bld = tidyr::replace_na(carbo_bld, median(carbo_bld, na.rm = TRUE)),
                                                         pot_bld = tidyr::replace_na(pot_bld, median(pot_bld, na.rm = TRUE)),
                                                         whg_bld = tidyr::replace_na(whg_bld, median(whg_bld, na.rm = TRUE)),
                                                         leg_bld = tidyr::replace_na(leg_bld, median(leg_bld, na.rm = TRUE)),
                                                         ref_bld = tidyr::replace_na(ref_bld, median(ref_bld, na.rm = TRUE)),
                                                         sugar_bld = tidyr::replace_na(sugar_bld, median(sugar_bld, na.rm = TRUE)),
                                                         sumref_bld = tidyr::replace_na(sumref_bld, median(sumref_bld, na.rm = TRUE)),
                                                         fruit_bld = tidyr::replace_na(fruit_bld, median(fruit_bld, na.rm = TRUE)),
                                                         veg_bld = tidyr::replace_na(veg_bld, median(veg_bld, na.rm = TRUE)),
                                                         act_bld = tidyr::replace_na(act_bld, median(act_bld, na.rm = TRUE)),
                                                         alco_bld = tidyr::replace_na(alco_bld, median(alco_bld, na.rm = TRUE)),
                                                         ahei_bld = tidyr::replace_na(ahei_bld, median(ahei_bld, na.rm = TRUE)),
                                                         aheinowgr_bld = tidyr::replace_na(aheinowgr_bld, median(aheinowgr_bld, na.rm = TRUE)),
                                                         calor_bld = tidyr::replace_na(calor_bld, median(calor_bld, na.rm = TRUE)),
                                                         alco_g=statar::xtile(alco_bld, n = 5),
                                                         act_g=statar::xtile(act_bld, n = 5),
                                                         ahei_g=statar::xtile(ahei_bld, n = 5),
                                                         aheinowgr_g=statar::xtile(aheinowgr_bld, n = 5),
                                                         prot_g=statar::xtile(prot_bld, n = 5),
                                                         calor_g=statar::xtile(calor_bld, n = 5)) %>% as.data.frame()

temp_cohort<-temp_cohort %>% mutate(dtdxdb2=as.numeric(dtdxdb2),
                                    dtdxcvd=as.numeric(dtdxcvd),
                                    dtdxchd=as.numeric(dtdxchd),
                                    dtdxstr=as.numeric(dtdxstr),
                                    dtdth=as.numeric(dtdth),
                                    cutoff=as.numeric(cutoff),
                                    blddate=as.numeric(blddate),
                                    dtdxdb2=case_when(is.na(dtdxdb2) ~ 9999, T ~ dtdxdb2),
                                    dtdth=case_when(is.na(dtdth) ~ 9999, T ~ dtdth),
                                    dtdxcvd=case_when(is.na(dtdxcvd) ~ 9999, T ~ dtdxcvd))  


# temp_cohort<-temp_cohort %>% mutate(db_time=case_when(db_out==1 ~ dtdxdb2-blddate,T ~ cutoff-blddate),
#                                    dth_time=case_when(dth_out==1 ~ dtdth-blddate,T ~ cutoff-blddate),
#                                    cvd_time=case_when(cvd_out==1 ~ dtdxcvd-blddate,T ~ cutoff-blddate),
#                                    chd_time=case_when(chd_out==1 ~ dtdxchd-blddate,T ~ cutoff-blddate),
#                                    str_time=case_when(str_out==1 ~ dtdxstr-blddate,T ~ cutoff-blddate))

temp_cohort<-temp_cohort %>% filter(!duplicated(newid)) # id no duplicate

table(temp_cohort$db2case,temp_cohort$cohort)
table(temp_cohort$cvdcase,temp_cohort$cohort)
# temp_cohort <- temp_cohort %>% mutate(across(c(carbo_bld, prot_bld, pot_bld, leg_bld, whg_bld, ref_bld),
#                                             ~ . * 4 / calor_bld, .names = "{.col}_energy"))


### Cox regression
fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld, data=temp_cohort)
fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld  +
            dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
            factor(act_g) + factor(alco_g)+ factor(hbp_bld) + 
            factor(chol_bld) , data=temp_cohort)

fit=coxph(Surv(tdb2,db2case) ~ y_predictSD + cohort + age_bld  + bmi_young +
            dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
            factor(act_g) + factor(alco_g)+ factor(hbp_bld) + 
            factor(chol_bld) , data=temp_cohort)


R<-data.frame(variable='model 1',HR=summary(fit)$conf.int[,1],LL=summary(fit)$conf.int[,3],
              UL=summary(fit)$conf.int[,4],P=summary(fit)$coefficients[,5])


R<-R %>% mutate(HRCI=paste0(round(HR,3),"(",round(LL,3),"-",round(UL,3),")"),Pvalue=round(P,3))


R
# total carb
#                   variable        HR        LL        UL            P               HRCI Pvalue
#y_predictSD       model 1 1.0861146 1.0326846 1.1423090 1.329368e-03 1.086(1.033-1.142)  0.001

### quintile HR 
#  y_predictSD quintile (1, 2, 3, 4, 5)
temp_cohort$y_predictSD_quintile <- ntile(temp_cohort$y_predictSD, 5)
temp_cohort$y_predictSD_quintile <- factor(temp_cohort$y_predictSD_quintile, levels = 1:5)
fit1_q <- coxph(Surv(tdb2, db2case) ~ y_predictSD_quintile + cohort + age_bld, data = temp_cohort)
fit2_q <- coxph(Surv(tdb2, db2case) ~ y_predictSD_quintile + cohort + age_bld + 
                  dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                  factor(act_g) + factor(alco_g) + factor(hbp_bld) + factor(chol_bld), 
                data = temp_cohort)
dim(summary(fit1_q)$conf.int)
dim(summary(fit1_q)$coefficients)
rownames(summary(fit1_q)$coefficients)

quintile_rows <- grep("y_predictSD_quintile", rownames(summary(fit1_q)$conf.int))


fit1_confint <- summary(fit1_q)$conf.int[quintile_rows, ]
fit1_coeff <- summary(fit1_q)$coefficients[quintile_rows, ]

R_fit1 <- data.frame(
  variable = paste0("Q", 2:5),  #  Q2-Q5
  HR = fit1_confint[, 1],       #  HR
  LL = fit1_confint[, 3],       # 
  UL = fit1_confint[, 4],       # 
  P = fit1_coeff[, 5]           # P 
)

R_fit1 <- R_fit1 %>%
  mutate(HRCI = paste0(round(HR, 2), "(", round(LL, 2), "-", round(UL, 2), ")"),
         Pvalue = round(P, 3))


quintile_rows_fit2 <- grep("y_predictSD_quintile", rownames(summary(fit2_q)$conf.int))
fit2_confint <- summary(fit2_q)$conf.int[quintile_rows_fit2, ]
fit2_coeff <- summary(fit2_q)$coefficients[quintile_rows_fit2, ]
R_fit2 <- data.frame(
  variable = paste0("Q", 2:5),  #  Q2-Q5 
  HR = fit2_confint[, 1],       #  HR
  LL = fit2_confint[, 3],        
  UL = fit2_confint[, 4],       
  P = fit2_coeff[, 5]           #  P 
)

R_fit2 <- R_fit2 %>%
  mutate(HRCI = paste0(round(HR, 2), "(", round(LL, 2), "-", round(UL, 2), ")"),
         Pvalue = round(P, 3))

R_combined <- rbind(R_fit1, R_fit2)
R_combined

# variable       HR        LL       UL            P            HRCI Pvalue
# y_predictSD_quintile2        Q2 1.351673 1.144260 1.596683 3.919153e-04  1.35(1.14-1.6)  0.000
# y_predictSD_quintile3        Q3 1.556979 1.323857 1.831151 8.792176e-08 1.56(1.32-1.83)  0.000
# y_predictSD_quintile4        Q4 1.617107 1.375894 1.900609 5.484793e-09  1.62(1.38-1.9)  0.000
# y_predictSD_quintile5        Q5 1.702375 1.450507 1.997977 7.383672e-11     1.7(1.45-2)  0.000
# y_predictSD_quintile21       Q2 1.241872 1.050407 1.468237 1.122540e-02 1.24(1.05-1.47)  0.011
# y_predictSD_quintile31       Q3 1.380696 1.172222 1.626246 1.121710e-04 1.38(1.17-1.63)  0.000
# y_predictSD_quintile41       Q4 1.351887 1.146296 1.594351 3.407801e-04 1.35(1.15-1.59)  0.000
# y_predictSD_quintile51       Q5 1.368453 1.160595 1.613539 1.901398e-04 1.37(1.16-1.61)  0.000


# P for trend
quintile_medians <- temp_cohort %>%
  group_by(y_predictSD_quintile) %>%
  summarise(median_value = median(y_predictSD, na.rm = TRUE))
temp_cohort <- temp_cohort %>%
  left_join(quintile_medians, by = c("y_predictSD_quintile" = "y_predictSD_quintile"))


fit1_trend <- coxph(Surv(tdb2, db2case) ~ median_value + cohort + age_bld, data = temp_cohort)
fit2_trend <- coxph(Surv(tdb2, db2case) ~ median_value + cohort + age_bld + 
                      dbfh + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                      factor(act_g) + factor(alco_g) + factor(hbp_bld) + factor(chol_bld), 
                    data = temp_cohort)
p_for_trend_fit1 <- summary(fit1_trend)$coefficients["median_value", "Pr(>|z|)"]
p_for_trend_fit2 <- summary(fit2_trend)$coefficients["median_value", "Pr(>|z|)"]
p_for_trend_fit1
p_for_trend_fit2

