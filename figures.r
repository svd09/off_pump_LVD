# figures for the OPCAB paper.
# 12/28/2020.

library(easypackages)
libraries(
    "tidyverse", "survival", "survminer",
    "rms", "Hmisc", "ggsci","lattice","ggthemes",
    "extrafont"
)

# get the data.

df <- read_csv("E:/opcab/opcab_paper/long_term_results/dataset_latest.csv")


# make the survival time to years.

df$surv_years <- (1 + df$survival_days) / 365.24

describe(df$surv_years)

# figure 1A. KM plot for survival for the whole group.

df$opcab <- factor(df$opcab)

s <- survfit(Surv(surv_years, died) ~ opcab, data = df)

survdiff(Surv(surv_years, died) ~ opcab, data = df)

# survival at time points.

summary(s, times = c(0, 2, 4, 6, 8, 10))

summary(s, times = c(0, 1,5, 10))

tiff("E:/opcab/opcab_paper/propensity_model/figures/figure1a.tiff",
    height = 5, width = 8, units = "in", res = 1200
)

ggsurvplot(s,
    data = df,
    risk.table = T,
    conf.int = T,
    break.x.by = 1,
    xlim = c(0, 10),
    palette = c("blue","red"),
    surv.scale = "percent",
    legend.labs = c("ONCAB", "OPCAB")
)

dev.off()

# now to do the same figure for the matched data.

# get the matched dataaset first.

df_m <- read_csv("E:\\opcab\\opcab_paper\\propensity_model\\prop_match_data.csv")

glimpse(df_m)

# now this is the propensity matched dataset.

df_m$surv_years <- (1 + df_m$survival_days) / 365.24

df_m$opcab <- factor(df_m$opcab)

s_m <- survfit(Surv(surv_years, died) ~ opcab, data = df_m)

summary(s_m, times = c(0, 2, 4, 8, 10))

summary(s_m, times = c(0,1,5,10))


survdiff(Surv(surv_years, died) ~ opcab, data = df_m)


tiff("E:/opcab/opcab_paper/propensity_model/figures/figure1b.tiff",
    height = 5, width = 8, units = "in", res = 1200
)

ggsurvplot(s_m,
    data = df,
    risk.table = T,
    conf.int = T,
    break.x.by = 1,
    xlim = c(0, 10),
    palette = c("blue","red"),
    surv.scale = "percent",
    legend.labs = c("ONCAB", "OPCAB")
)

dev.off()

# figure for the propensity scores model.
# mirror histogram.

# for that we need to first again calculate the PS scores.


formula <- opcab ~ age_at_surgery + gender + bmi + diabetes +
    dialysis_dependence + systemic_hypertension + smoker + hyperlipidemia +
    copd + peripheral_vascular_disease + preop_lvef + elective + urgent + emergency +
    left_main_disease + preop_mi + preop_cva + prior_cardiac_surgery + logistic_euroscore +
    coronary_disease_detail

psmodel <- glm(formula, data = df, family = "binomial"(link = "logit"))

df$ps <- fitted.values(psmodel)

# now we can separate the ps scores according to category.

opcab <- df %>% filter(opcab == 1)
oncab <- df %>% filter(opcab == 0)

# now to create the mirror histogram for the data.

library(ggsci)

mypal <- pal_jama(palette = c("default"), alpha = 1)

tiff("E:/opcab/opcab_paper/propensity_model/figures/ps_score.tiff",
    height = 5, width = 7, units = "in", res = 1200
)

par(mfrow = c(2, 1))
par(family = 'Ariel')

# Make the plot
par(mar = c(0, 5, 3, 3))
hist(opcab$ps,
    main = "",
    xlim = c(0, 1), ylab = "OPCAB",
    xlab = "", ylim = c(0, 30), xaxt = "n", las = 1,
    col = "red",
    breaks = 100
)
par(mar = c(5, 5, 0, 3))
hist(oncab$ps,
    main = "",
    xlim = c(0, 1), ylab = "ONCAB",
    xlab = "Propensity Score distribution",
    ylim = c(30, 0), las = 1, col = "blue",
    breaks = 100
)

dev.off()

# create the love plot for the matched data.

# obtain the data from the cobalt library and then
# use the lattice to plot the data.
# plot like from the MAG/SAG paper can be done

# get the data_dc for dotchart.

d = read_csv("E:/opcab/opcab_paper/propensity_model/data_dc.csv")

d2 = d %>% pivot_longer(!Variable,
                  names_to = "type")

d2$Variable = factor(d2$Variable, ordered = T)

a = ggplot(data = d2, aes(x = value, y = Variable, group = type,
                      color = type)) + 
  
  geom_point() + 
  geom_vline(xintercept = 0.1, col = "red", linetype = 2) + 
  geom_vline(xintercept = -0.1, col = "red", linetype = 2) +
  ylab(" ") + xlab("Standardised Difference value") + theme_clean() + scale_color_aaas() +
  geom_vline(xintercept = 0, col = "blue", linetype = 2)

a2 = a + theme_clean() + theme(legend.position = "none") + 
  scale_x_continuous(breaks = c(-0.2,-0.1,0,0.1,0.2,0.3,0.4)) +
    theme(text=element_text(family="Arial", size = 12)) +
  scale_color_manual(values = c("black","gray"))

a2

tiff(  filename = "E:/opcab/opcab_paper/propensity_model/figures/lp.tiff",
      
       height = 10, width = 6,
       unit = "in", res = 1200)

a2

dev.off()


