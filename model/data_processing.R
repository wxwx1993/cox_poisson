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
setwd("/Users/apple/Dropbox/Documents/GitHub/cox_poisson")
file.list <- list.files("synthetic_medicare/merged_synthetic_mortality", pattern = ".csv", full.names = TRUE)

var.list <- c("year", "fipscounty", "BENE_SEX_IDENT_CD", "BENE_RACE_CD", "age", "state", "region", 
              "dead", "pm25", "mean_bmi", "smoke_rate", "blk_pct", "median_household_income", 
              "median_house_value", "poverty", "no_grad", "population_density", "owner_occupied",
              "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", # more 
              "DESYNPUF_ID", "BENE_BIRTH_DT", "BENE_DEATH_DT")

rename.list <- c("sex" = "BENE_SEX_IDENT_CD", 
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
#mean(df$entry_age_break >= 65) # 0.837
age_break <- c(seq(min(df$entry_age_break), 100, by=5), max(df$entry_age_break))
age_label <- c("25-29", "30-34", "35-39", "40-44", "45-49", 
               "50-54", "55-59", "60-64", "65-69", "70-74", 
               "75-79", "80-84", "85-89", "90-94", "95-99", ">= 100")
df$entry_age_break <- cut(df$entry_age_break, breaks = age_break, right = FALSE, labels = age_label)

### follow_up_year is the gap between yeat_at_entry and calendar year
df$followup_year <- df$year - df$entry_year
df$followup_year_plus_one <- df$followup_year + 1

###################################################################################################
# Generate count data for each individual characteristics and follow-up year
df$time_count <- df$followup_year_plus_one - df$followup_year # always 1
dead_personyear <- aggregate(cbind(df$dead, df$time_count), 
                             by=list(df$fipscounty,
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
                         by=list(df$fipscounty,
                                 df$year,
                                 df$sex,
                                 df$race,
                                 df$entry_age_break,
                                 df$followup_year),
                         FUN = min) 


aggregate_data <- merge(dead_personyear,
                        confounders,
                        by = c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5", "Group.6"))
names(aggregate_data)[1:8] <- c("fipscounty", "year", "sex", "race", "entry_age_break", "followup_year", "dead", "time_count")
aggregate_data <- subset(aggregate_data[complete.cases(aggregate_data), ])

### quick check 
sum(aggregate_data$dead) == sum(df$dead) 
sum(aggregate_data$time_count) == dim(df)[1]

save(aggregate_data, file = "data/aggregate_data.RData")
