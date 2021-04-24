# code for propensity matched model.
# propensity matching 1:1 for OPCAB and ONCAB.
# use code from the lita/bita paper for ease.

# get the libraries needed.

library(pacman)
p_load(
    tidyverse, survival, rms, Hmisc, survminer,
    survRM2, relsurv
)


library(easypackages)
libraries(c(
    "tidyverse", "haven", "readxl",
    "survival", "lubridate", "WeightIt",
    "cobalt", "rms", "Hmisc", "survminer",
    "cmprsk", "splines"))

library(cobalt)
library(Matching)
library(MatchIt)
library(WeightIt)
library(twang)
library(MatchThem)
library(naniar)
library("MatchIt")
library("lmtest")
library("sandwich")
library("boot")
library("survival")

# get the dataset here.

df <- read_csv("D:/opcab/opcab_paper/long_term_results/dataset_latest.csv")

glimpse(df)

#THE ORIGINAL PAPER WAS SENT WITHOUT ADDING THE SURGEON VOLUME DATA.
#DF2 <- DF WILL GET THE ORIGINAL RESULTS.
#REMOVE TYPE FROM THE PROPENSITY MODEL.
# to add surgeon volume here...

dfs <- read_excel("D:/opcab/opcab_paper/propensity_model/raw_data.xlsx",sheet = 1)

dfs2 <- dfs %>% dplyr::select(Pat_ID, Surgeon_name)


dfs2 <- dfs2[!duplicated(dfs2$Pat_ID), ]


s_count <- dfs2 %>% group_by(Surgeon_name) %>% summarise(total = n())

# list of low:- 

s_l <- s_count %>% filter(total <= 30)

s_m <- s_count %>% filter(total > 30 & total <= 60)

s_h <- s_count %>% filter(total > 60)

# now add the surgeon_names to the main dataset.

names(dfs2) <- tolower(names(dfs2))

dfs2$pat_id = as.numeric(dfs2$pat_id)

df2 <- left_join(df, dfs2, by = "pat_id")


df2$type <- with(df2, ifelse(surgeon_name %in% s_l$Surgeon_name, 1,
                             ifelse(surgeon_name %in% s_m$Surgeon_name, 2, 3)))


df2 %>% count(type)

# now to see missing information for the variables to be used for the propensity score model.

vars <- c("age_at_surgery",
    "gender", "bmi",
    "diabetes", "dialysis_dependence", "systemic_hypertension", "smoker",
    "hyperlipidemia", "copd", "peripheral_vascular_disease", "preop_lvef",
    "elective", "urgent", "emergency", "left_main_disease", "preop_mi",
    "preop_cva", "prior_cardiac_surgery", "logistic_euroscore", "opcab", "type")

miss_var_summary(df2[, c(vars)])

summary(df2[, c(vars)])

# while none missing, preop_mi has -9 and bmi needs to change.

df2$preop_mi[df2$preop_mi == -9] <- 001

df2$bmi[df2$bmi == 94] <- 49

df2%>% count(coronary_disease_detail)

# now to use all these variables and create a model for propensity score.


formula <- opcab ~ age_at_surgery + gender + bmi + diabetes +
    dialysis_dependence + systemic_hypertension + smoker + hyperlipidemia +
    copd + peripheral_vascular_disease + preop_lvef + elective + type +  
  urgent + emergency + left_main_disease + preop_mi + preop_cva + prior_cardiac_surgery + logistic_euroscore +
    coronary_disease_detail

psmodel <- glm(formula, data = df2, family = "binomial"(link = "logit"))

df2$ps <- fitted.values(psmodel)

# using the MatchIt package to do the matching and then use the cobalt library
# to do the balance checks.


library(mosaic)

favstats(ps ~ opcab, data = df2)


# opcab        min        Q1    median        Q3       max      mean         sd   n
# 1     0 0.08759807 0.2982455 0.3508001 0.4266387 0.7205048 0.3640819 0.09584981 719
# 2     1 0.17017071 0.3334291 0.3963809 0.4743310 0.7580399 0.4077491 0.10461494 442


results = broom::tidy(psmodel)

write_csv(results,
    "D:/opcab/opcab_paper/propensity_model/table_model.csv")

gm <- matchit(data = df2, formula = formula,
    method = "nearest", caliper = 0.2,
    replace = F, ratio = 1)


love.plot(gm)

cobalt::bal.tab(gm, un = T)

lp <- cobalt::bal.tab(gm, un = T)

data = lp$Balance %>% tbl_df()

Variable = c("Distance","PS", "Sex","BMI","DM","ESRD",
             "HTN","Smoking","Hyperlipidemia","COPD",
             "PAD","LVEF","Elective","Urgent","Emergent",
             "LMCA stenosis","Prior MI","Prior Stroke",
             "Prior OHS","EUROSCORE", "CAD")

data$Variable = Variable

data2 = data %>% dplyr::select(Variable, Diff.Un, Diff.Adj)

write_csv(data2,
          "D:/opcab/opcab_paper/propensity_model/data_dc.csv")

    # the balance looks good.
    # now to extract the data-set using 1:1 matching.
    # get the data-set.


df_m <- match.data(gm)

df_m %>% count(opcab)

# now to save this dataaset for further analysis.

write_csv(df_m,
    "G:/opcab/opcab_paper/propensity_model/prop_match_data.csv")

# now to see the cox model and then km plot for the matched data.

df_m$surv_years <- (df_m$survival_days + 1) / 365.24

km <- survfit(Surv(surv_years, died) ~ opcab, data = df_m)


ggsurvplot(km,
    censor.size = 0,
    conf.int = T,
    risk.table = T,
    xlim = c(0, 10),
    surv.scale = "percent",
    legend.label = c("ONCAB", "OPCAB"))


# get the Cox model with robust standard errors.

coxm <- coxph(Surv(surv_years, died) ~ opcab + incomplete_revascularisation, data = df_m)


cox_st <- coxph(Surv(surv_years, died) ~ opcab + strata(subclass) + incomplete_revascularisation,
data = df_m)

summary(cox_st)

summary(coxm) # this is the result after matching 1:1.


# doing a sub-group analysis for OPCAB/ONCAB and incomplete grafting.
# for OPCAB and for ONCAB separately.

# using the adjusted data.


tableone::CreateCatTable(vars = "opcab",
        data = df_m,
        strata = c("incomplete_revascularisation"))

# for incomp == 1.


d1 = df_m %>% filter(incomplete_revascularisation == 1)


m1 = coxph(Surv(surv_years, died) ~ opcab, data = d1)

summary(m1)

# now adjusted for all the other variables along with opcab.

formula <- Surv(surv_years, died) ~ age_at_surgery + gender + bmi + diabetes +
dialysis_dependence + systemic_hypertension + smoker + hyperlipidemia +
copd + peripheral_vascular_disease + preop_lvef + elective + urgent + emergency +
left_main_disease + preop_mi + preop_cva + prior_cardiac_surgery + logistic_euroscore +
 opcab 


m1_adj = coxph(formula, data = d1)


d2 = df_m %>% filter(incomplete_revascularisation == 0)

m2 = coxph(Surv(surv_years, died) ~ opcab, data = d2)

summary(m2)

m2_adj = coxph(formula, data = d2)

summary(m2_adj)


# create new term for interaction.


df_m$inc_opcab = df_m$incomplete_revascularisation*df_m$opcab


formula <- Surv(surv_years, died) ~ age_at_surgery + gender + bmi + diabetes +
  dialysis_dependence + systemic_hypertension + smoker + hyperlipidemia +
  copd + peripheral_vascular_disease + preop_lvef + elective + urgent + emergency +
  left_main_disease + preop_mi + preop_cva + prior_cardiac_surgery + logistic_euroscore +
  opcab + incomplete_revascularisation + inc_opcab


m_interact = coxph(formula, data = df_m)

# some more code to fill the tables.

df <- read_csv("D:/opcab/opcab_paper/propensity_model/prop_match_data.csv")


# male gender.

tableone::CreateCatTable(vars = "gender",
                         data = df,
                         strata = c("opcab"))

# LOS median IQR and p-value:

h <- tableone::CreateContTable(vars = "hospitallenthof_stay",
                          data = df,
                          strata = c("opcab"),
                          testNonNormal = T)

print(h, nonnormal = c(" hospitallenthof_stay"))

library(mosaic)

favstats(hospitallenthof_stay ~ opcab, data = df)



# survival years and pt-years follow-up:

df$surv_years <- df$survival_days/365.24

mean(df$surv_years) * 1161

df %>% count(died)

319/(mean(df$surv_years) * 1161)