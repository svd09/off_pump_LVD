library(tidyverse);library(haven);
library(readxl)

df = 
read_excel("G:\\opcab\\opcab_paper\\long_term_results\\final_data.xlsx",
sheet = 1)

df2 = df[!duplicated(df$Pat_ID), ]



df3 = df2 %>% select(Pat_ID,
                     Survival_days,
                     Died)

names(df3) = tolower(names(df3))

df4 = df3 %>% rename(
  surv_d_n = survival_days,
  died_n = died
)


l = 
read_csv("E:\\opcab\\opcab_paper\\long_term_results\\dataset_latest.csv")

glimpse(l)


df4$pat_id = as.numeric(df4$pat_id)


l2 = left_join(l, df4, by = "pat_id")

write_csv(l2, 
"E:\\opcab\\opcab_paper\\long_term_results\\l2_lt_latest.csv" )


table(l2$died_n)