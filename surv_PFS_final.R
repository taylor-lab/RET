#Code for generating survival plots for figures 2b and extended data 2a, 2b, and 7b

library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)
library(survival)

require(survminer)

#read in data
df_all <- data.frame(fread("~/Documents/17-256/new_analysis/PFS_table.csv", header = TRUE, sep = ",", fill = TRUE))

#Shedding - Figure 2B
df_Shedding <- df_all[df_all$Shedding != "REMOVE",]
df_Shedding <- df_Shedding[df_Shedding$Shedding != "GERMLINE",]

fitShedding =survfit(Surv(Weeks, Progression_Yes_No) ~ Shedding_Binary, df_Shedding)
surv_median(fitShedding)
pN = ggsurvplot(fitShedding, xlim = c(0, 120), legend.labs = c("Shedding=No", "Shedding=Yes"))$plot
pN + labs(x="Weeks",
          y="Progression Free Survival")
summary(coxph(Surv(Weeks, Progression_Yes_No) ~ Shedding_Binary, data = df_Shedding))

#TP53 status - Extended Data Figure 2A
df_TP53 <- df_all[df_all$TP53 != "REMOVE",]
fitTP53 =survfit(Surv(Weeks, Progression_Yes_No) ~ TP53_Binary, df_TP53)
surv_median(fitTP53)
pN = ggsurvplot(fitTP53, xlim = c(0, 120), legend.labs = c("TP53=No", "TP53=Yes"))$plot
pN + labs(x="Weeks",
          y="Progression Free Survival")
summary(coxph(Surv(Weeks, Progression_Yes_No) ~ TP53_Binary, data = df_TP53))

#prev MKI Anal - Extended Data Figure 2B
df_MKI <- df_all[df_all$Treatment_Naive != "REMOVE",]
fitMKI =survfit(Surv(Weeks, Progression_Yes_No) ~ Previous_MKI_Binary, df_MKI)
surv_median(fitMKI)
pN = ggsurvplot(fitMKI, xlim = c(0, 120), legend.labs = c("Previous MKI=No", "Previous MKI=Yes"))$plot
pN + labs(x="Weeks",
          y="Progression Free Survival")
summary(coxph(Surv(Weeks, Progression_Yes_No) ~ Previous_MKI_Binary, data = df_MKI))

#PI3K status - Extended Data Figure 7B
df_PI3K <- df_all[df_all$PI3K != "REMOVE",]
fitPI3K =survfit(Surv(Weeks, Progression_Yes_No) ~ PI3K_Binary, df_PI3K)
surv_median(fitPI3K)
pN = ggsurvplot(fitPI3K, xlim = c(0, 120), legend.labs = c("PI3K=No", "PI3K=Yes"))$plot
pN + labs(x="Weeks",
          y="Progression Free Survival")
summary(coxph(Surv(Weeks, Progression_Yes_No) ~ PI3K_Binary, data = df_PI3K))
