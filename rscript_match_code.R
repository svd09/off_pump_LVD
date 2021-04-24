# further analysis using the matched dataset.

library(easypackages)
libraries(
    "tidyverse", "tableone",
    "rms", "Hmisc", "rstpm2", "naniar",
    "survminer", "lintr", "ResourceSelection"
)





df_m <-
    read_csv("D:/opcab/opcab_paper/propensity_model/prop_match_data.csv")

# now to look at the intraoperative details in the matched dataset.

# table 1 with matched dataset.

vars <- c(
    "age_at_surgery",
    "gender", "bmi",
    "diabetes", "dialysis_dependence", "systemic_hypertension", "smoker",
    "hyperlipidemia", "copd", "peripheral_vascular_disease", "preop_lvef",
    "elective", "urgent", "emergency", "left_main_disease", "preop_mi",
    "preop_cva", "prior_cardiac_surgery", "coronary_disease_detail",
    "additive_euroscore", "logistic_euroscore", "opcab"
)


factors <- c(
    "gender",
    "diabetes", "dialysis_dependence", "systemic_hypertension", "smoker",
    "hyperlipidemia", "copd", "peripheral_vascular_disease",
    "elective", "urgent", "emergency", "left_main_disease", "preop_mi",
    "preop_cva", "prior_cardiac_surgery",
    "opcab", "coronary_disease_detail"
)



t1_matched <- tableone::CreateTableOne(
    vars = vars,
    factorVars = factors,
    data = df_m,
    strata = c("opcab")
)


t1_match <- print(t1_matched,
    nonnormal = c(
        "age_at_surgery", "bmi",
        "logistic_euroscore", "additive_euroscore"
    )
)

write.csv(
    t1_match,
    "G:/opcab/opcab_paper/propensity_model/t1_match.csv"
)


# from table 1, it appears that the baseline data is now very well matched.
# in the earlier code, forgot CAD.

df_m %>% count(coronary_disease_detail)

t_c = tableone::CreateCatTable(vars = "coronary_disease_detail",
                               data = df_m,
                               strata = c("opcab"))


t_c

# look at operative details for the matched patients.

vars <- c(
    "opcab", "length_of_surgery",
    "incomplete_revascularisation",
    "no_of_revascularizedvessels",
    "bima", "radialis", "veins"
)

factors <- c(
    "opcab",
    "incomplete_revascularisation",

    "bima", "radialis", "veins"
)

t2 <- tableone::CreateTableOne(
    vars = vars,
    data = df_m,
    strata = c("opcab"),
    factorVars = factors
)

t2

t2_match <- print(t2,
    nonnormal = c("length_of_surgery")
)


# now to save this and then work on the variables that
# are important for only OPCAB or ONCAB.

write.csv(
    t2_match,
    "G:/opcab/opcab_paper/propensity_model/t2_match.csv"
)

# now to look at the outcomes for the matched data.

# conversion to onpump.

opcab <- df_m %>% filter(opcab == 1)

tableone::CreateCatTable(
    vars = "conversion_to_onpump",
    data = opcab
)

# Overall
# n                            430
# conversion_to_onpump = 1 (%) 20 (4.7)

# for the oncab patients now.

oncab <- df_m %>% filter(opcab == 0)

tableone::CreateCatTable(
    vars = c(
        "onpump_beating_heart",
        "on_pump_cross_clamp"
    ),
    data = oncab
)


# Overall
# n                            430
# onpump_beating_heart = 1 (%) 117 (27.2)
# on_pump_cross_clamp = 1 (%)  312 (72.6)

time <- tableone::CreateContTable(vars = c("bypass_time",
                                   "cross_clamp_time"),
                          data = df_m[df_m$opcab == 0, ])

print(time, nonnormal = c("bypass_time","cross_clamp_time"))

# now to look at the postoperative outcomes for the matched data.

vars <- c(
    "postop_low_cardiac_output",
    "postop_iabp", "postop_ecmo", "postop_mi", "postop_resuscitation",
    "postop_cardiac_arrhythmia", "postop_redo_heart", "postop_reexploration_bleeding",
    "postop_symptomatic_transitory_ps", "pstop_cva", "postop_sepsis",
    "postop_new_dialysis", "ekrbcs", "tkplatelets", "ffp", "postop_respiratory_failure",
    "re_intubation", "tracheotomy", "hospitallenthof_stay",
    "died_30days"
)

factorvars <- c(
    "postop_low_cardiac_output",
    "postop_iabp", "postop_ecmo", "postop_mi", "postop_resuscitation",
    "postop_cardiac_arrhythmia", "postop_redo_heart", "postop_reexploration_bleeding",
    "postop_symptomatic_transitory_ps", "pstop_cva", "postop_sepsis",
    "postop_new_dialysis", "postop_respiratory_failure",
    "re_intubation", "tracheotomy",
    "died_30days"
)

t3_outcome_matched <- tableone::CreateTableOne(
    vars = vars,
    factorVars = factorvars,
    data = df_m, strata = c("opcab")
)

t3_outcome_matched <- print(t3_outcome_matched,
    nonnormal = c("hospitallenthof_stay")
)

write.csv(
    t3_outcome_matched,
    "G:/opcab/opcab_paper/propensity_model/t3_matched_outcome.csv"
)

# glm model for 30 day mortality adjusting for incomplete repeat_revascularization
# and number of grafts.

m1 <- glm(died_30days ~ opcab + incomplete_revascularisation,
    data = df_m, family = "binomial"(link = "logit")
)

summary(m1)

exp(coef(m1))
exp(confint(m1))

# function to obtain robust estimates and SE for a glm model:

# this function can be used for obtaining robust estimates for GLM models.


robust_est <- function(fitglm) {
    require(sandwich)
    cov.m1 <- vcovHC(fitglm, type = "HC0")

    std.err <- sqrt(diag(cov.m1))

    q.val <- qnorm(0.975)

    r.est <- cbind(
        Estimate = exp(coef(fitglm)),
        "Robust SE" = std.err,
        z = (coef(fitglm) / std.err),
        "Pr(>|z|) " = 2 * pnorm(abs(coef(fitglm) / std.err), lower.tail = FALSE),
        LL = exp(coef(fitglm) - q.val * std.err),
        UL = exp(coef(fitglm) + q.val * std.err)
    )
    print(r.est)
}

robust_est(m1)

hoslem.test(df_m$died_30days, fitted(m1), g = 10)

# p = 0.99 for hosmer lemeshow test shows good fit for glm model


# now to look at the long term outcome for the matched data


df_m$surv_years <- (df_m$survival_days + 1) / 365.24

surv_m <- survfit(Surv(surv_years, died) ~ opcab, data = df_m)

ggsurvplot(surv_m,
    risk.table = T, censor.size = 0,
    conf.int = T, xlim = c(0, 8),
    break.x.by = 2
)


# now to do the model for late survival.


cox_matched <- coxph(Surv(surv_years, died) ~ opcab +
    incomplete_revascularisation, data = df_m)


summary(cox_matched)

cox.zph(cox_matched)


# create a newdataset for 4 patients and then plot the fitted model.


newdata <- tibble(
    opcab = c(0, 0, 1, 1),
    incomplete_revascularisation = c(0, 1, 0, 1)
)

new_predict <- survfit(cox_matched,
    newdata = newdata,
    conf.int = T, conf.type = "log-log"
)

glimpse(new_predict)

newp2 <- broom::tidy(new_predict)

glimpse(newp2)

plot(new_predict,
    xlim = c(0, 8),

    col = c("red", "red", "blue", "blue"),
    lty = c(1, 1, 2, 2)
)




# opcab vs oncab + beating vs oncab + xclamp.
# rstmp2 model.
# after looking at this, we are not going to include this in the paper.

df_m2 <- df_m %>% filter(opcab == 0)


surv_m2 <- survfit(Surv(surv_years, died) ~
onpump_beating_heart, data = df_m2)

ggsurvplot(surv_m2,
    risk.table = T, censor.size = 0,
    conf.int = T, xlim = c(0, 8),
    break.x.by = 2
)


cox1 <- coxph(Surv(surv_years, died) ~
onpump_beating_heart, data = df_m2)



vars <- c(
    "onpump_beating_heart",
    "age_at_surgery",
    "gender", "bmi",
    "diabetes", "dialysis_dependence", "systemic_hypertension", "smoker",
    "hyperlipidemia", "copd", "peripheral_vascular_disease", "preop_lvef",
    "elective", "urgent", "emergency", "left_main_disease", "preop_mi",
    "preop_cva", "prior_cardiac_surgery",
    "additive_euroscore", "logistic_euroscore", "opcab"
)


factors <- c(
    "gender",
    "diabetes", "dialysis_dependence", "systemic_hypertension", "smoker",
    "hyperlipidemia", "copd", "peripheral_vascular_disease",
    "elective", "urgent", "emergency", "left_main_disease", "preop_mi",
    "preop_cva", "prior_cardiac_surgery",
    "opcab"
)



t1_matched2 <- tableone::CreateTableOne(
    vars = vars,
    factorVars = factors,
    data = df_m2,
    strata = c("onpump_beating_heart")
)

t1_matched2


# see how incomplete grafting changes the km curve.

inc <- survfit(Surv(surv_years, died) ~ incomplete_revascularisation,
    data = df_m
)

ggsurvplot(inc,
    risk.table = T, conf.int = T,
    censor.size = 0,
    xlim = c(0, 8),
    break.x.by = 2
)


survdiff(Surv(surv_years, died) ~ incomplete_revascularisation,
    data = df_m
)

