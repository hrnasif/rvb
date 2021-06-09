## HERS Data Clean up
# Data sourced from the Heart and Estrogen/Progestin Study (HERS, Hulley et al., 1998)
# Downloaded from Table 11.8 at:
# https://regression.ucsf.edu/second-edition/data-examples-and-problems
hers_raw <- read.csv("../hersdata.csv")

exclusions <- which(
  is.na(hers_raw$sbp) | # No systolic blood pressure data
    is.na(hers_raw$bmi) | # No BMI data
    is.na(hers_raw$htnmeds) | # No data on whether medication is taken
    is.na(hers_raw$age) | # No data on age
    is.na(hers_raw$visit)) # No data on visit number

hers_excl <- hers_raw[-exclusions,]

id <- hers_excl$pptid
unique_id <- unique(id) # 2031 unique patients

# Define variables
response <- as.numeric(hers_excl$sbp > 140) # Response
htn <- hers_excl$htnmeds # Already defined as 0/1 indicator
bmi <- rvb::standardize(hers_excl$bmi) # standardize BMI
age <- rvb::standardize(hers_excl$age) # standardize age
visit <- hers_excl$visit*2/5 - 1 # Visit number coded

rvb_hers <- data.frame(id = id,
                       response = response,
                       htn = htn,
                       bmi = bmi,
                       age = age,
                       visit = visit)

usethis::use_data(rvb_hers, overwrite = TRUE)
