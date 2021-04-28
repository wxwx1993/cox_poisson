# library("yaml")
# The following objects are masked from 'package:data.table':
# hour, isoweek, mday, minute, month, quarter, second, wday, week, yday, year
# library(parallel)
# library("foreign")

library(dplyr)
library(readr)
library(lubridate)
library(survival)
library(gnm)
###################################################################################################
setwd("/Users/mac/Documents/GitHub/survival")
file.list <- list.files("data/merged_synthetic_mortality", pattern = ".csv", full.names = TRUE)

var.list <- c("year", "fipscounty", "BENE_SEX_IDENT_CD", "BENE_RACE_CD", "age", "state", "region", 
              "dead", "pm25", "mean_bmi", "smoke_rate", "blk_pct", "median_household_income", 
              "median_house_value", "poverty", "no_grad", "population_density", "owner_occupied",
              "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", # more 
              "DESYNPUF_ID", "BENE_BIRTH_DT", "BENE_DEATH_DT")

rename.list <- c("zip" = "fipscounty", "sex" = "BENE_SEX_IDENT_CD", 
                 "race" = "BENE_RACE_CD", "statecode" = "state", 
                 "pm25_ensemble" = "pm25", "pct_blk" = "blk_pct", 
                 "medhouseholdincome" = "median_household_income", 
                 "medianhousevalue" = "median_house_value", "education" = "no_grad",
                 "popdensity" = "population_density", "pct_owner_occ" = "owner_occupied",
                 "id" = "DESYNPUF_ID", "birthdate" = "BENE_BIRTH_DT", "deathdate" = "BENE_DEATH_DT")

df <- file.list %>% lapply(read_csv) %>% lapply(subset, select = var.list) %>% bind_rows %>% rename(rename.list)

### formatting 
df <- df[complete.cases(df[, setdiff(names(df), c("id", "birthdate", "deathdate"))]), ]
df$birthdate <- ymd(df$birthdate)
df$deathdate <- ymd(df$deathdate)
df$id <- as.character(df$id)

### for each individual, for the first year in record
entry_year <- df %>% group_by(id) %>% summarise(entry_year = min(year))
df <- merge(df, entry_year, on="id", all.x = T)
rm(entry_year)

### year_at_entry = birth_year + year_at_entry
df$entry_age_break <- df$entry_year - year(df$birthdate)
age_break <- c(seq(25, 95, by=5), 105)
age_label <- c("25-29", "30-34", "35-39", "40-44", "45-49", 
               "50-54", "55-59", "60-64", "65-69", "70-74", 
               "75-79", "80-84", "85-89", "90-94", "95-101")
df$entry_age_break <- cut(df$entry_age_break, breaks = age_break, right = FALSE, labels = age_label)

### follow_up_year is the gap between yeat_at_entry and calendar year
df$followup_year <- df$year - df$entry_year
df$followup_year_plus_one <- df$follow_up_year + 1

### quick check 
summary(df$year - year(df$birthdate) - df$age)
summary(year(df$deathdate) - year(df$birthdate) - df$age) 
print("age column is calendar-year age (death age as well). should no use.") 
summary(as.factor(df$followup_year))
print("about 92% id retained in 2010 since 2008")
# 0       1       2 
# 2211358 2176107 2033946 
summary(as.factor(df$entry_year))
print("few new entry in 2009 and 2010")
# 2008    2009    2010 
# 6419904    1344     163 

summary(df$entry_age_break)
print("dominated by elder >= 65")

save(df, file = paste0("data/df.RData"))
# load("data/df.RData")





###################################################################################################
# Generate count data for each individual characteristics and follow-up year
df$time_count <- df$followup_year_plus_one - df$followup_year # always 1??? 
dead_personyear <- aggregate(cbind(df$dead, df$time_count), 
                             by=list(df$zip,
                                     df$year,
                                     df$sex,
                                     df$race,
                                     df$entry_age_break,
                                     df$followup_year),
                             FUN = sum)

confounder.list <- c("pm25_ensemble", "mean_bmi", "smoke_rate", "pct_blk",
                     "medhouseholdincome", "medianhousevalue", "poverty", 
                     "education", "popdensity", "pct_owner_occ", "summer_tmmx", 
                     "winter_tmmx", "summer_rmax", "winter_rmax", "region")

confounders <- aggregate(df[,confounder.list],
                         by=list(df$zip,
                                 df$year,
                                 df$sex,
                                 df$race,
                                 df$entry_age_break,
                                 df$followup_year),
                         FUN = min) 
# df[,confounder.list] %>% group_by(zip, ...) %>% summarise(=[1])


aggregate_data <- merge(dead_personyear,
                        confounders,
                        by = c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5", "Group.6"))
names(aggregate_data)[1:8] <- c("zip", "year", "sex", "race", "entry_age_break", "followup_year", "dead", "time_count")
aggregate_data <- subset(aggregate_data[complete.cases(aggregate_data), ])

### quick check 
sum(aggregate_data$dead) == sum(df$dead) 
sum(aggregate_data$time_count) == dim(df)[1]

save(aggregate_data, file = "data/aggregate_data.RData")

### ????
table(as.factor(df$year), as.factor(df$entry_year))

###################################################################################################
# Cox Proportional Hazard
# ds4 <- df[sample(nrow(df), 600000*2), ]

start_time <- Sys.time()
Cox_raw <- coxph(Surv(followup_year, followup_year_plus_one, dead) ~ 
                   pm25_ensemble + mean_bmi + smoke_rate + pct_blk +
                   medhouseholdincome + medianhousevalue +
                   poverty + education + popdensity + pct_owner_occ +
                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                   as.factor(year) +
                   as.factor(region) +
                   strata(as.factor(entry_age_break)) +
                   strata(as.factor(sex)) +
                   strata(as.factor(race)),
                 data = df, ties = c("efron"), na.action = na.omit)
end_time <- Sys.time()
print(end_time - start_time)
Cox <- summary(Cox_raw)
save(Cox, file = paste0("data/Cox.RData"))



###################################################################################################
# Cox-equvalent conditional Poisson Regression
gnm_raw <- gnm(dead ~ pm25_ensemble + 
                 mean_bmi + smoke_rate + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ +
                 summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region) +
                 offset(log(time_count)),
               eliminate = (as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(followup_year)),
               data = aggregate_data,
               family = poisson(link = "log"))
Poisson <- summary(gnm_raw)
# exp(10 * Poisson$coefficients[1])
save(Poisson, file = "data/Poisson.RData")

# we allow multiple row for each indivisual? for? 
# read_yaml("/Users/mac/Documents/GitHub/survival/data/merged_synthetic_mortality/schema.yml")
# d08 = read.csv("/Users/mac/Documents/GitHub/survival/data/merged_synthetic_mortality/synthetic_with_confounders_2008.csv")
