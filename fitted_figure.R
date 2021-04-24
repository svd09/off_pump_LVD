# fitted Cox model plot for OPCAB/ONCAB and incomp revasc.

# get the libraries and the matched dataset.



library(easypackages)
libraries(c(
  "tidyverse", "haven", "readxl",
  "survival",  "rms", "Hmisc", "survminer","ggsci"
))

mycol = pal_jama(palette = c("default"), alpha = 1)(4)

df_m <- read_csv("E:/opcab/opcab_paper/propensity_model/prop_match_data.csv")

# create cox model for opcab and incomp revasc.


df_m$surv_years <- (df_m$survival_days + 1) / 365.24

coxm <- coxph(Surv(surv_years, died) ~ opcab + 
                incomplete_revascularisation, data = df_m)



newdata = tibble(opcab = c(1,1,0,0),
                 incomplete_revascularisation = c(0,1,0,1)
                 )

fm <- survfit(coxm, newdata = newdata)




tiff("E:/opcab/opcab_paper/propensity_model/fit_model_new.tiff",
    
    height = 5,
    width = 7,
    units = "in",
  res = 1200)

plot(fm,xlim = c(0,10), col = c("red","red","blue","blue"), 
     lty = c(1,2,1,2),
     ylab = "Predicted Survival", 
     xlab = "Time:Years", frame.plot = F)

dev.off()

# to obtain survival estimates at 5 and 10 years.


newdata_5 = tibble(opcab = c(1,1,0,0),
                 incomplete_revascularisation = c(0,1,0,1),
                 surv_years = c(5,5,5,5),
                 died = c(0,0,0,0))

pm_5 <- predict(coxm, newdata = newdata_5,
                type = "expected")

newdata_10 = tibble(opcab = c(1,1,0,0),
                   incomplete_revascularisation = c(0,1,0,1),
                   surv_years = c(10,10,10,10),
                   died = c(1,1,1,1))


pm_10 <- predict(coxm, newdata = newdata_10,
                type = "expected")


