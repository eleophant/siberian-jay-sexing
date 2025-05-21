####
#
# siberian jay sexing
#
####

# packages ----
library(tidyverse)

# data input ----
sexing = read_csv("siberian_jay_sexing.csv")
metadata_1 = read_csv("df_metadata_sib.csv")
metadata_2 = read_csv("fasteris_samples_metadata.csv")

# sexing data ----
sexing |>  select(sample_name) |> head(n = 3) #last 8 characters are ring number

# make a new column with ring numbers 
sexing = sexing |>
  mutate(ring_number = str_sub(sample_name, -7))

# assign sex (if allele1 = allele2, then male)
sexing = sexing |> 
  mutate(sex = ifelse(allele1 == allele2, 'm', 'f'))

# _ tidy ----
# identify ring_number duplicates and discrepancies
sexing |> 
  group_by(ring_number) |> 
  filter(n() > 1) |>   #keep duplicated ring_numbers
  summarise(
    sex_values = n_distinct(sex),
    sex_consistent = n_distinct(sex) == 1, #there are 34 duplicates
    .groups = "drop") |> #ungroup
  filter(sex_consistent == FALSE) #9 rows are discrepant

# 9/25 duplicates are discrepant: 36% error rate

# my metadata ----
# get all ring numbers from metadata 1 and 2
metadata_1 |> select(ring_number)
metadata_2 |> select(ring_number)

# new df with ring_number and sex
metadata_1$ring_number = as.character(metadata_1$ring_number) 

metadata_1_sexing = metadata_1 |> mutate(sex_original = sex) |> 
  select(ring_number, sex_original)

metadata_1_sexing = metadata_1_sexing |> 
  left_join(sexing, by = "ring_number") |> 
  mutate(sex_may2025 = sex) |> 
  select(ring_number, sex_original, sex_may2025)

metadata_1_sexing |> 
  filter(sex_original != sex_may2025) #no discrepancies

metadata_1_sexing |> 
  filter(!is.na(sex_may2025)) #only 49 of the 91 have sexing data

# repeat with metadata_2
# new df with ring_number and sex
metadata_2$ring_number = as.character(metadata_2$ring_number) 

metadata_2_sexing = metadata_2 |> 
  left_join(sexing, by = "ring_number") |> 
  mutate(sex_may2025 = sex) |> 
  select(ring_number, sex_may2025)

metadata_2_sexing |> 
  filter(!is.na(sex_may2025)) #only 33 of the 71 have sexing data


