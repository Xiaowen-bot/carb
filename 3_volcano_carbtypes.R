
#-------------------------------------------------------------------------------
#  December 7, 2024
#  The purpose of this program is to:
#   (1) Individual regression: metabolites~diet (refined grains, vegetables, fruits, potatoes, legumes)
#   (2) draw volcano plots
#  Figure 3
#-------------------------------------------------------------------------------

names(temp_lvs)
met_list=as.character(colnames(temp_lvs[,2:294]))
results_lvs=data.frame(met=met_list)
met_aligned=results_lvs %>% pull(met)

for (i in 1:293) {
  fit=lm(as.formula(paste0(met_aligned[i] ,"~ refined_grain_r + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_rgrain"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_rgrain"]=summ$coefficients[2,2]
  results_lvs[i, "p_rgrain"]=summ$coefficients[2,4]
  results_lvs[i, "n_rgrain"]=nrow(model.frame(fit))
  fit=lm(as.formula(paste0(met_aligned[i], "~ refined_grain_r + cohort + ageyr + white + factor(bmi_group) + factor(neversmoking) + 
  factor(ahei_g) + factor(calor_g) + factor(act_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_rgrain_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_rgrain_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_rgrain_adj"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i],"~ vegetable_r + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_vegetable"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_vegetable"]=summ$coefficients[2,2]
  results_lvs[i, "p_vegetable"]=summ$coefficients[2,4]
  results_lvs[i, "n_vegetable"]=nrow(model.frame(fit))
  fit=lm(as.formula(paste0(met_aligned[i], "~ vegetable_r + cohort + ageyr  + white + factor(bmi_group) + factor(neversmoking) +
  factor(alco_g)+ factor(ahei_noveg_g) + factor(calor_g) + factor(act_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_vegetable_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_vegetable_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_vegetable_adj"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i],"~ fruit_r + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_fruit"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_fruit"]=summ$coefficients[2,2]
  results_lvs[i, "p_fruit"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ fruit_r + cohort  + factor(bmi_group) + ageyr  + white + factor(neversmoking) + factor(ahei_nofrt_g) + factor(calor_g) + 
                           factor(act_g) + factor(alco_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_fruit_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_fruit_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_fruit_adj"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i],"~ legume_r + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_legume"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_legume"]=summ$coefficients[2,2]
  results_lvs[i, "p_legume"]=summ$coefficients[2,4]
  results_lvs[i, "n_legume"]=nrow(model.frame(fit))
  fit=lm(as.formula(paste0(met_aligned[i], "~ legume_r + cohort + ageyr  + white + factor(bmi_group) + factor(neversmoking) +
  factor(alco_g)+ factor(ahei_noleg_g) + factor(calor_g) + factor(act_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_legume_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_legume_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_legume_adj"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i],"~ potato_r + cohort")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_potato"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_potato"]=summ$coefficients[2,2]
  results_lvs[i, "p_potato"]=summ$coefficients[2,4]
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ potato_r + cohort  + factor(bmi_group) + ageyr  + white + factor(neversmoking) + factor(ahei_g) + factor(calor_g) + 
                           factor(act_g) + factor(alco_g)")), data=temp_lvs)
  summ=summary(fit)
  results_lvs[i, "beta_potato_adj"]=round(summ$coefficients[2,1],4)
  results_lvs[i, "se_potato_adj"]=summ$coefficients[2,2]
  results_lvs[i, "p_potato_adj"]=summ$coefficients[2,4]
  
}

results_lvs<- results_lvs %>% mutate(
  rgrain_intake=case_when(beta_rgrain_adj>0 & p_rgrain_adj<0.05/293  ~ 'positive',
                          beta_rgrain_adj<0 & p_rgrain_adj<0.05/293  ~ 'negative',
                          T ~ 'insignificant'),
  vegetable_intake=case_when(beta_vegetable_adj>0 & p_vegetable_adj<0.05/293  ~ 'positive',
                            beta_vegetable_adj<0 & p_vegetable_adj<0.05/293  ~ 'negative',
                            T ~ 'insignificant'),
  fruit_intake=case_when(beta_fruit_adj>0 & p_fruit_adj<0.05/293  ~ 'positive',
                         beta_fruit_adj<0 & p_fruit_adj<0.05/293  ~ 'negative',
                         T ~ 'insignificant'),
  legume_intake=case_when(beta_legume_adj>0 & p_legume_adj<0.05/293  ~ 'positive',
                         beta_legume_adj<0 & p_legume_adj<0.05/293  ~ 'negative',
                         T ~ 'insignificant'),
  potato_intake=case_when(beta_potato_adj>0 & p_potato_adj<0.05/293  ~ 'positive',
                          beta_potato_adj<0 & p_potato_adj<0.05/293  ~ 'negative',
                          T ~ 'insignificant'))

f_lvs_match<-f_lvs %>% mutate(met=rownames(f_lvs)) %>% 
  dplyr::select("met",'method',"metabolite_name",'biochemical_name',"class_metabolon","sub_class_metabolon","super_class_metabolon")
results_lvs<-merge(results_lvs, f_lvs_match, by = 'met', all = FALSE, sort = TRUE)

write.csv(results_lvs, "/udd/n2xwa/carb_met/results_lvs_20250103.csv")


#--------------------------------volcano plot ----------------------------------  

results_vol <- results_lvs
alpha_adjusted <- 0.05 / 293

# refined grain

results_vol <- transform(results_vol, color = ifelse(rgrain_intake == "positive", "#E41A1C",
                                                     ifelse(rgrain_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(rgrain_intake == "positive" | rgrain_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_rgrain_adj < alpha_adjusted )

r <- ggplot(results_vol, 
            aes(x = beta_rgrain_adj, y = -log10(p_rgrain_adj), 
                color = color, label = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha_adjusted), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") + 
  geom_text_repel(data = significant_metabolites, aes(x = beta_rgrain_adj, y = -log10(p_rgrain_adj), label = metabolite_name), 
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
  labs(title = "Refined Grain Consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "Negative"))
r
ggsave("/udd/n2xwa/carb_met/volcano_rgrain.png", plot = r, width = 8, height = 8, dpi = 1200)


# select top ten: fruit
results_vol <- results_lvs
results_vol <- transform(results_vol, color = ifelse(fruit_intake == "positive", "#E41A1C",
                                                     ifelse(fruit_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(fruit_intake == "positive" | fruit_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_fruit_adj < alpha_adjusted )

p2 = ggplot(results_vol, aes(x = beta_fruit_adj, y = -log10(p_fruit_adj), 
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
  labs(title = "Fruit consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "negative"))  

# select top10
Up <- filter(results_vol, fruit_intake == 'positive') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_fruit_adj))
Down <- filter(results_vol, fruit_intake == 'negative') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_fruit_adj))

nudge_x_up = 0.005 - Up$beta_fruit_adj 
nudge_x_down = -0.005 - Down$beta_fruit_adj

p3 <- p2 + 
  geom_point(data = Up,aes(x = beta_fruit_adj, y = -log10(p_fruit_adj)),
             color = '#EB4232', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = beta_fruit_adj, y = -log10(p_fruit_adj), label = metabolite_name),
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
f <- p3 + 
  geom_point(data = Down,aes(x = beta_fruit_adj, y = -log10(p_fruit_adj)),
             color = '#2DB2EB', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Down,aes(x = beta_fruit_adj, y = -log10(p_fruit_adj), label = metabolite_name),
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



f
ggsave("/udd/n2xwa/carb_met/volcano_fruit.png", plot = f, width = 8, height = 8, dpi = 1200)


# select top ten: vegetables

# select top ten
results_vol <- results_lvs
results_vol <- transform(results_vol, color = ifelse(vegetable_intake == "positive", "#E41A1C",
                                                     ifelse(vegetable_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(vegetable_intake == "positive" | vegetable_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_vegetable_adj < alpha_adjusted )

p2 = ggplot(results_vol, aes(x = beta_vegetable_adj, y = -log10(p_vegetable_adj), 
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
  labs(title = "Vegetable consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "negative"))  

# select top10
Up <- filter(results_vol, vegetable_intake == 'positive') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_vegetable_adj))
Down <- filter(results_vol, vegetable_intake == 'negative') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_vegetable_adj))

nudge_x_up = 0.005 - Up$beta_vegetable_adj 
nudge_x_down = -0.005 - Down$beta_vegetable_adj

p3 <- p2 + 
  geom_point(data = Up,aes(x = beta_vegetable_adj, y = -log10(p_vegetable_adj)),
             color = '#EB4232', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = beta_vegetable_adj, y = -log10(p_vegetable_adj), label = metabolite_name),
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
v <- p3 + 
  geom_point(data = Down,aes(x = beta_vegetable_adj, y = -log10(p_vegetable_adj)),
             color = '#2DB2EB', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Down,aes(x = beta_vegetable_adj, y = -log10(p_vegetable_adj), label = metabolite_name),
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

v

ggsave("/udd/n2xwa/carb_met/volcano_vegetable.png", plot = v, width = 8, height = 8, dpi = 1200)

# legume
results_vol <- results_lvs
results_vol <- transform(results_vol, color = ifelse(legume_intake == "positive", "#E41A1C",
                                                     ifelse(legume_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(legume_intake == "positive" | legume_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_legume_adj < alpha_adjusted )

l <- ggplot(results_vol, 
            aes(x = beta_legume_adj, y = -log10(p_legume_adj), 
                color = color, label = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha_adjusted), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") + 
  geom_text_repel(data = significant_metabolites, aes(x = beta_legume_adj, y = -log10(p_legume_adj), label = metabolite_name), 
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
  labs(title = "Legume Consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "Negative"))
l
ggsave("/udd/n2xwa/carb_met/volcano_legume.png", plot = l, width = 8, height = 8, dpi = 1200)

# potato
results_vol <- results_lvs
results_vol <- transform(results_vol, color = ifelse(potato_intake == "positive", "#E41A1C",
                                                     ifelse(potato_intake == "negative", "#377EB8", "#d8d8d8")),
                         significance = ifelse(potato_intake == "positive" | potato_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_potato_adj < alpha_adjusted )

p <- ggplot(results_vol, 
            aes(x = beta_potato_adj, y = -log10(p_potato_adj), 
                color = color, label = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha_adjusted), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") + 
  geom_text_repel(data = significant_metabolites, aes(x = beta_potato_adj, y = -log10(p_potato_adj), label = metabolite_name), 
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
  labs(title = "Potato Consumption",
       x = "Adjusted-beta",
       y = "-log10(p-value)",
       color = "Metabolites") +
  scale_color_manual(values = c("#d8d8d8" = "#d8d8d8", "#E41A1C" = "#E41A1C", "#377EB8" = "#377EB8"),
                     labels = c("#d8d8d8" = "Insignificance", "#E41A1C" = "Positive", "#377EB8" = "Negative"))
p
ggsave("/udd/n2xwa/carb_met/volcano_potato.png", plot = p, width = 8, height = 8, dpi = 1200)

# Sub-arrangement 1: 3 plots in a single row
subplot1 <- ggarrange(p4, w, r, 
                      ncol = 3, 
                      labels = c("A", "B", "C"),
                      widths = c(0.35, 0.35, 0.3),  
                      common.legend = TRUE,
                      legend = "top")

# Sub-arrangement 2: 4 plots in a single row
f2 <- f + theme(legend.position = "none")
v2 <- v + theme(legend.position = "none")
l2 <- l + theme(legend.position = "none")
p2 <- p + theme(legend.position = "none")
subplot2 <- ggarrange(f2, v2, 
                      ncol = 2, 
                      labels = c("D", "E"),
                      widths = c(0.5, 0.5))
subplot3 <- ggarrange(l2, p2,
                      ncol = 2, 
                      labels = c("F", "G"),
                      widths = c(0.5, 0.5))
f_f <- ggarrange(subplot1, subplot2, subplot3,
                 nrow = 3,
                 heights = c(0.35, 0.35, 0.3))
ggsave(filename = "/udd/n2xwa/carb_met/figure_3.pdf", plot = f_f, 
       width = 12,#
       height = 16, #
       # units = "in",
       dpi = 800)

