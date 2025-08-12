####
#
#
# siberian_jays_analysis.R
#
#
####

# packages ----
library(tidyverse)
library(conflicted)
library(phyloseq)
library(fantaxtic)
library(data.table)
library(ggordiplots)
library(ggvegan)
library(eulerr)
library(microbiome)
library(viridis)
library(vegan)
library(compositions)
library(ALDEx2)
theme_set(theme_classic())
conflicted::conflicts_prefer(dplyr::filter())
conflicted::conflicts_prefer(dplyr::select())

# data input ----
df_metadata_sib = read_csv("df_metadata_sib_no_duplicates_filtered.csv")
df_otus_sib = read_csv("df_otus_sib_no_duplicates.csv")
taxonomy_table = read_csv("jay_otu_taxonomy_table.csv")
territory_sizes_f22 = read_csv("territories_f22.csv")

# check if any split samples (samples from the same bird on the same day)
df_metadata_sib |> distinct(ring_number, date, .keep_all = TRUE) #there are none

# generate otu matrix
df_otus_sib = df_otus_sib |> column_to_rownames("rowname")

# otu cleanup ----
# pull Chloroplast and Mitochondria
chlo_mito_seqs = taxonomy_table |> 
  filter(Order == "Chloroplast" | Family == "Mitochondria") |> 
  pull(...1)

# remove Chloroplast and Mitochondria from df_otus_sib
chlo_mito_seqs_valid = intersect(chlo_mito_seqs, colnames(df_otus_sib))

df_otus_sib_filtered = df_otus_sib |> 
  as.data.frame() |> 
  select(-all_of(chlo_mito_seqs_valid))

# check for dataset-wide singletons
asv_total_counts = colSums(df_otus_sib_filtered) #there are 2

# get dataset-wide singletons
singleton_asvs = asv_total_counts[asv_total_counts == 1]

# view singletons
names(singleton_asvs)

# remove singletons from df_otus
df_otus_sib_filtered_no_singletons = df_otus_sib_filtered |> 
  select(-all_of(singleton_asvs))

#write_csv(df_otus_sib_filtered_no_singletons, "df_otus_sib_filtered_no_singletons.csv")

# overwrite df_otus_sib for downstream analysis
df_otus_sib = df_otus_sib_filtered_no_singletons
rm(df_otus_sib_filtered, df_otus_sib_filtered_no_singletons)

# subset fall22 ----
# subset metadata
df_metadata_sib_f22 = df_metadata_sib |> filter(season == "fall2022")
df_metadata_sib_f22 = df_metadata_sib_f22 |> left_join(territory_sizes_f22, by = join_by(territory))

# subset otus
samples_f22 = df_metadata_sib |> filter(season == "fall2022") |> pull(study_id)
df_otus_sib_f22 = df_otus_sib |> as.data.frame() |> rownames_to_column() |> filter(rowname %in% samples_f22) |> column_to_rownames("rowname")
df_otus_sib_f22 = df_otus_sib_f22[,-(which(colSums(df_otus_sib_f22) == 0))]

# phyloseq ----
otu = df_otus_sib |> t()
otu = otu_table(otu, taxa_are_rows = TRUE)
tax = taxonomy_table |> column_to_rownames("...1")
tax = tax_table(as.matrix(tax))
metadata = df_metadata_sib |> column_to_rownames("study_id")
metadata = sample_data(metadata)
siberian = phyloseq(otu, tax, metadata)

# sample sizes ----
df_metadata_sib |> group_by(season) |> summarise(n = n()) #winter2022: 12, summer2022: 8, fall2022: 49, winter2023: 13

df_metadata_sib |> group_by(area) |> summarise(n = n()) #man 34, reivo 48

territorial_samples = df_metadata_sib |> group_by(territory) |> summarise(n = n()) |> mutate(obs = n) 
territorial_samples |> summarise(mean(obs), max(obs), min(obs)) #mean = 1.95, max = 7, min = 1

df_metadata_sib |> group_by(sample_type) |> summarise(n = n()) #cloacal 53, faecal 29

df_metadata_sib |> group_by(ring_number) |> summarise(n = n()) #73 individuals

df_metadata_sib |> group_by(age) |> summarise(n = n()) #adult 47, juvenile 35

df_metadata_sib |> group_by(breeding_status) |> summarise(n = n()) #breeder 31, non-breeder 50, NA 1

# reads & ASVs ----
# total number of reads
df_otus_sib |> sum() #3,748,491

# number of reads per sample
mean(rowSums(df_otus_sib)) #45,713.3
sd(rowSums(df_otus_sib)) #62,898.55
max(rowSums(df_otus_sib)) #401,095
min(rowSums(df_otus_sib)) #37

# number of ASVs
df_otus_sib |> ncol() #3331

# number of ASVs per sample
# number of non-zero rows for each column is the number of ASVs for that sample
# samples need to be in columns
total_asvs = df_otus_sib |> 
  t() |> as.data.frame() |> 
  summarise(across(starts_with("KCG"), ~ sum(.x != 0))) |> 
  pivot_longer(cols = everything()) 

mean(total_asvs$value) #87.87805
sd(total_asvs$value) #78.35191
max(total_asvs$value) #405
min(total_asvs$value) #8

# relative abundances ---- 
# normalize number of reads using median sequencing depth
total = median(sample_sums(siberian))
standf = function(x, t = total) round(t * (x / sum(x)))
a_sib_norm = transform_sample_counts(siberian, standf)

# plot relative abundances
plot_bar(a_sib_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") #+ facet_grid(~ area)

# dominant taxa ----
fantaxtic::top_taxa(siberian, n = 5, tax_level = "Phylum")
fantaxtic::top_taxa(siberian, n = 5, tax_level = "Class")
fantaxtic::top_taxa(siberian, n = 5, tax_level = "Order")
fantaxtic::top_taxa(siberian, n = 5, tax_level = "Family")
fantaxtic::top_taxa(siberian, n = 5, tax_level = "Genus")
fantaxtic::top_taxa(siberian, n = 5, tax_level = "Species")

# plot dominant taxa
top_asv = fantaxtic::top_taxa(siberian, n_taxa = 10)
plot_nested_bar(ps_obj = top_asv$ps_obj,
                top_level = "Phylum",
                nested_level = "Genus")

# alpha diversity ----

# if df_metadata_sib_no_duplicates_filtered.csv was loaded, these steps have already been done

# _ adiv data ----
adiv = plot_richness(siberian, measures = c("Observed", "Chao1", "Shannon", "Simpson"))

# store alpha diversity measures as new variables
alphadiv = data.table(adiv$data)

# create tibble with adiv measures for each sample
alphadiv = alphadiv |> 
  select(mcmaster_sample_id, samples, variable, value) |> 
  pivot_wider(names_from = variable, values_from = value) |> 
  mutate(mcmaster_sample_id = as.character(mcmaster_sample_id))

# bind to existing metadata
alphadiv$mcmaster_sample_id = as.numeric(alphadiv$mcmaster_sample_id)

# merge to metadata
df_metadata_sib = df_metadata_sib |> 
  left_join(alphadiv, by = "mcmaster_sample_id")

write_csv(df_metadata_sib2, "df_metadata_sib_no_duplicates_filtered.csv")

# _ time ----
df_metadata_sib |>
  ggplot(aes(x = julian, y = Observed, color = territory))+
  geom_point(size = 1) +
  ggtitle("Alpha diversity (Observed ASVs) over time")

kruskal.test(Observed ~ julian, data = df_metadata_sib) #Kruskal-Wallis chi-squared = 36.258, df = 38, p-value = 0.5502

# _ season ----
season_order <- c('winter2022', 'summer2022', 'fall2022','winter2023') 

df_metadata_sib$season <- factor(df_metadata_sib$season, levels = season_order)

df_metadata_sib |>
  ggplot(aes(x = season, y = Observed, color = season)) +
  geom_point() +
  geom_boxplot(size = 1)
# +ggtitle("Alpha diversity (Observed ASVs) per season")

df_metadata_sib |>
  ggplot(aes(x = season, y = Shannon, color = season)) +
  geom_boxplot(size = 1)

kruskal.test(Observed ~ season, data = df_metadata_sib) #Kruskal-Wallis chi-squared = 9.7141, df = 3, p-value = 0.02116

kruskal.test(Shannon ~ season, data = df_metadata_sib) #Kruskal-Wallis chi-squared = 4.4473, df = 3, p-value = 0.217

pairwise.wilcox.test(df_metadata_sib$Observed, df_metadata_sib$season, p.adjust.method = "bonferroni") #only diffs bw winter22 & fall22, winter22 & winter23; no diffs when testing Shannon

# _ area ----
df_metadata_sib |>
  ggplot(aes(x = area, y = Observed, color = area))+
  geom_boxplot(size = 1) #+
#ggtitle("Alpha diversity (Observed ASVs) per area")

df_metadata_sib |>
  ggplot(aes(x = area, y = Shannon, color = area))+
  geom_boxplot(size = 1) #+
#ggtitle("Alpha diversity (Shannon PD) per area")

# test area
kruskal.test(Observed ~ area, data = df_metadata_sib) #Kruskal-Wallis chi-squared = 1.0145, df = 1, p-value = 0.3138

kruskal.test(Shannon ~ area, data = df_metadata_sib) #Kruskal-Wallis chi-squared = 0.42178, df = 1, p-value = 0.5161

# _ breeding status ----
df_metadata_sib |>
  filter(breeding_status != "NA") |> 
  ggplot(aes(x = breeding_status, y = Observed, color = age))+
  geom_boxplot(size = 1) +
  #ggtitle("Alpha diversity (Observed ASVs) vs breeding status") +
  xlab("Breeding status") + ylab("Observed ASVs")

df_metadata_sib |>
  filter(breeding_status != "NA") |> 
  ggplot(aes(x = breeding_status, y = Shannon, color = age))+
  geom_boxplot(size = 1) +
  #ggtitle("Alpha diversity (Shannon PD) \nvs breeding status (Siberian jays)") +
  xlab("Breeding status")

# test alpha diversity by breeding status
kruskal.test(Observed ~ breeding_status, data = df_metadata_sib) #Kruskal-Wallis chi-squared = 1.4755, df = 1, p-value = 0.2245

kruskal.test(Shannon ~ breeding_status, data = df_metadata_sib) #Kruskal-Wallis chi-squared = 0.084972, df = 1, p-value = 0.7707

# _ age ----
df_metadata_sib |>
  ggplot(aes(x = age, y = Observed, colour = age))+
  geom_boxplot(size = 1) +
  xlab("Age") + ylab("Observed ASVs") #+
#ggtitle("Alpha diversity (Observed ASVs) vs age") +

kruskal.test(Shannon ~ age, data = df_metadata_sib) #Shannon KW chi-squared = 0.36564, df = 1, p = 0.5454

kruskal.test(Observed ~ age, data = df_metadata_sib) #Observed KW chi-squared = 0.60555, df = 1, p = 0.4365

# _ territory f22 ----
df_metadata_sib |>
  filter(season == "fall2022") |> 
  ggplot(aes(x = territory, y = Observed, color = territory))+
  geom_point(size = 1) +
  ggtitle("Alpha diversity (Observed ASVs) per territory \n(Fall 2022)")

kruskal.test(Shannon ~ territory, data = df_metadata_sib_f22) #p = 0.57 & 0.39 for Obs & Shannon

# _ size f22 ----
df_metadata_sib_f22 |>
  ggplot(aes(x = group_size, y = Observed, color = julian))+
  geom_point(size = 1) +
  ggtitle("Alpha diversity (Observed ASVs) by group size \n(Fall 2022)")

kruskal.test(Shannon ~ group_size, data = df_metadata_sib_f22) #p = 0.56 & 0.56 for Obs & Shannon

# beta diversity ----
# NMDS ----
# ordinate using NMDS on robust Aitchison distances
nmds_siberian = metaMDS(df_otus_sib, distance = "robust.aitchison", autotransform = FALSE) #autotransform needs to be false, otherwise it applies the Wisconsin transformation, which is not necessary as we are already using the robust Aitchison transformation

# add metadata to ordination
nmds_fort_siberian = fortify(nmds_siberian) |> 
  filter(score == "sites") |>
  rename(study_id = label) |> 
  left_join(df_metadata_sib)

# plot ordination
ggplot(nmds_fort_siberian, aes(x = NMDS1, y = NMDS2, colour = season)) + geom_point() + geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") #+
  theme(legend.position = "none")

# spiderplot
gg_ordiplot(nmds_siberian, groups = df_metadata_sib$season, hull = FALSE, label = TRUE, spiders = TRUE,
            ellipse = FALSE, pt.size = 2, plot = TRUE)

# _ fall22 ----
# ordinate using NMDS on robust Aitchison distances
nmds_siberian_f22 = metaMDS(df_otus_sib_f22, distance = "robust.aitchison", autotransform = FALSE) #autotransform needs to be false, otherwise it applies the Wisconsin transformation, which is not necessary as we are already using the robust Aitchison transformation

# add metadata to ordination
nmds_fort_siberian_f22 = fortify(nmds_siberian_f22) |> 
  filter(score == "sites") |>
  rename(study_id = label) |> 
  left_join(df_metadata_sib_f22)

# plot ordination
ggplot(nmds_fort_siberian_f22, aes(x = NMDS1, y = NMDS2, colour = area)) + geom_point() + geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") + theme_classic() + theme(legend.position = "none")

# db-RDA ----
# NB: order of rows in df_metadata & df_otus needs to match
df_otus_sib = df_otus_sib |> column_to_rownames("rowname")

# _ time pcnm ----
# rda on julian
mod_day = rda(df_otus_sib ~ julian, data = df_metadata_sib)
anova(mod_day) #p = 0.188
RsquareAdj(mod_day) #R^2 = 0.00601124

# generate pcnm
pcnm_time = df_metadata_sib |> select(julian) |> compositions::dist(method = "euclidean") |> pcnm()
df_metadata_sib_time = df_metadata_sib |> cbind(pcnm_time[["vectors"]])

# rda on pcnm
mod_time = rda(df_otus_sib ~ PCNM1 + PCNM2 + PCNM3 + PCNM4, data = df_metadata_sib_time)
anova(mod_time) #p = 0.06491308
RsquareAdj(mod_time) #R^2 = 0.01633714

# _ season ---- 
rda_sib_season = capscale(formula = df_otus_sib ~ season, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_season) #F = 2.4539, p = 0.001
RsquareAdj(rda_sib_season) #adj R^2 = 0.04018907

rda_scores_sib_season_df = rda_sib_season |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

ggplot(rda_scores_sib_season_df, aes(x = CAP1, y = CAP2, colour = season)) + geom_point(size = 2) + geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted")

gg_ordiplot(rda_sib_season, groups = df_metadata_sib$season, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 5.74%, CAP2 1.68%

# _ dietary season ----
rda_sib_diet_season = capscale(formula = df_otus_sib ~ dietary_season, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_diet_season) #p = 0.001
RsquareAdj(rda_sib_diet_season) #adj R^2 = 0.03577286

rda_scores_sib_diet_season_df = rda_sib_diet_season |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

ggplot(rda_scores_sib_diet_season_df, aes(x = CAP1, y = CAP2, colour = dietary_season, shape = sample_type)) +
  geom_point()

# _ sample type ----
rda_sib_type = capscale(formula = df_otus_sib ~ sample_type, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_type) #p = 0.017
RsquareAdj(rda_sib_type) #R^2 = 0.005495911

rda_scores_sib_type_df = rda_sib_type |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |>
  mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

ggplot(rda_scores_sib_type_df, aes(x = CAP1, y = MDS1, colour = season, shape = sample_type)) + geom_point() #sample type seems to be confounded by season

# _ space pcnm ----
# compute euclidean distances from xy coordinates, generate pcnm, then run RDA
# make latitude and longitude numeric
df_metadata_sib$X = as.numeric(df_metadata_sib$X)
df_metadata_sib$Y = as.numeric(df_metadata_sib$Y)

# compute pcnm
pcnm_dist = df_metadata_sib |> select(X, Y) |> dist(method = "euclidean") |> pcnm()
df_metadata_sib_pcnm = df_metadata_sib |> cbind(pcnm_dist[["vectors"]])

mod_spatial = rda(df_otus_sib ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13, data = df_metadata_sib_pcnm)
anova(mod_spatial) #p = 0.912

# _ location Y ----
rda_sib_y = capscale(formula = df_otus_sib ~ Y, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_y) #p = 0.255

# _ area ----
rda_sib_area = capscale(formula = df_otus_sib ~ area, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude) 
anova(rda_sib_area) #p = 0.449

rda_scores_sib_area_df = rda_sib_area |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

ggplot(rda_scores_sib_area_df, aes(x = CAP1, y = MDS1, colour = area, shape = season)) + geom_point(size = 2) #area is confounded by season

# _ territory ----
rda_sib_territory = capscale(formula = df_otus_sib ~ territory, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude) 
anova(rda_sib_territory) #0.062 - 0.077
RsquareAdj(rda_sib_territory) #adj R^2 = 0.07607335

rda_scores_sib_territory_df = rda_sib_territory |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

ggplot(rda_scores_sib_territory_df, aes(x = CAP1, y = CAP2, colour = territory)) + geom_point(size = 2)

# _ ring ----
df_metadata_sib$ring_number = as.character(df_metadata_sib$ring_number)
rda_sib_ring = capscale(formula = df_otus_sib ~ ring_number, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_ring) #p = 0.659

# _ age ----
rda_sib_age = capscale(formula = df_otus_sib ~ age, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_age) #p = 0.008
RsquareAdj(rda_sib_age) #adj R^2 = 0.005761799

rda_scores_sib_age_df = rda_sib_age |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

ggplot(rda_scores_sib_age_df, aes(x = CAP1, y = MDS1, colour = season, shape = age)) + geom_point() #age is confounded by season

gg_ordiplot(rda_sib_age, groups = df_metadata_sib$age, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 = 1.8% (no CAP2)

# _ age fine ----
rda_sib_age_fine = capscale(formula = df_otus_sib ~ age_fine, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_age_fine) #p = 0.157

# _ breeding status ----
## all breeders are adults, but non-breeders are both
rda_sib_breeding = capscale(formula = df_otus_sib ~ breeding_status, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_breeding) #p = 0.001
RsquareAdj(rda_sib_breeding) #adj R^2 = 0.01128518

rda_scores_sib_breeding_df = rda_sib_breeding |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

ggplot(rda_scores_sib_breeding_df, aes(x = CAP1, y = MDS1, colour = season, shape = breeding_status)) + geom_point() #breeding status is confounded by season

# _ model selection ----
# run null and full models
mod0_sib = capscale(df_otus_sib ~ 1, data = df_metadata_sib, distance = "robust.aitchison", na.action = na.exclude)

mod1_sib = capscale(formula = df_otus_sib ~ season + breeding_status + age, data = df_metadata_sib, distance = "robust.aitchison", na.action = na.exclude) #omitted time PCNM 
anova(mod1_sib) #p = 0.001
RsquareAdj(mod1_sib) #ajd R^2 = 0.03627631

# model selection on R^2 and p-values
step_r2_sib = ordiR2step(mod0_sib, scope = formula(mod1_sib), perm.max = 999, na.action = na.exclude)
#season is the only significant predictor: adj R^2 = 0.051095955

# fall 2022 ----

# _ time pcnm ----
# generate temporal pcnm
pcnm_time_f22 = df_metadata_sib_f22 |> select(julian) |> dist(method = "euclidean") |> pcnm()
df_metadata_sib_time_f22 = df_metadata_sib_f22 |> cbind(pcnm_time_f22[["vectors"]])

# rda temporal pcnm
mod_time_f22 = rda(df_otus_sib_f22 ~ PCNM1 + PCNM2 + PCNM3 + PCNM4, data = df_metadata_sib_time_f22)
anova(mod_time_f22) #p = 0.199

# _ sample type ----
rda_sib_type_f22 = capscale(formula = df_otus_sib_f22 ~ sample_type, data = df_metadata_sib_f22,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_type_f22) #p = 0.247

# _ space pcnm ----
pcnm_dist_f22 = df_metadata_sib_f22 |> select(X, Y) |> dist(method = "euclidean") |> pcnm()
df_metadata_sib_pcnm_f22 = df_metadata_sib_f22 |> cbind(pcnm_dist_f22[["vectors"]])

mod_spatial_f22 = rda(df_otus_sib_f22 ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13, data = df_metadata_sib_pcnm_f22)
anova(mod_spatial_f22) #p = 0.769

# _ location Y ----
rda_sib_y_f22 = capscale(formula = df_otus_sib_f22 ~ Y, data = df_metadata_sib_f22,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_y_f22) #p = 0.392

# _ area ----
rda_sib_area_f22 = capscale(formula = df_otus_sib_f22 ~ area, data = df_metadata_sib_f22,  distance = "robust.aitchison", na.action = na.exclude) 
anova(rda_sib_area_f22) #0.191

# _ territory ----
rda_sib_territory_f22 = capscale(formula = df_otus_sib_f22 ~ territory, data = df_metadata_sib_f22,  distance = "robust.aitchison", na.action = na.exclude) 
anova(rda_sib_territory_f22) #F = 1.3473, p = 0.023
RsquareAdj(rda_sib_territory_f22) #adj R^2 = 0.1684623

rda_scores_sib_territory_df_f22 = rda_sib_territory_f22 |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib_f22, by = "study_id")

ggplot(rda_scores_sib_territory_df_f22, aes(x = CAP1, y = CAP2, colour = territory, shape = age)) + geom_point(size = 2)

df_metadata_sib_f22 |> group_by(territory, ring_number) |> summarise(n = n()) |> print(n = 48) #1 individual was sampled twice

# repeat this on 48 unique individuals
df_metadata_sib_f22_uniq = df_metadata_sib_f22 |> distinct(ring_number, .keep_all = TRUE)

samples_f22_uniq = df_metadata_sib_f22_uniq |> pull(rowname)

df_otus_sib_f22_uniq = df_otus_sib_f22 |> rownames_to_column() |> filter(rowname %in% samples_f22_uniq) |> column_to_rownames("rowname")

rda_sib_territory_f22_uniq = capscale(formula = df_otus_sib_f22_uniq ~ territory, data = df_metadata_sib_f22_uniq,  distance = "robust.aitchison", na.action = na.exclude) 
anova(rda_sib_territory_f22_uniq) #F = 1.3685, p = 0.03
RsquareAdj(rda_sib_territory_f22_uniq) #adj R^2 = 0.1800064

rda_scores_sib_territory_df_f22_uniq = rda_sib_territory_f22_uniq |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib_f22_uniq, by = "study_id")

# plots
ggplot(rda_scores_sib_territory_df_f22_uniq, aes(x = CAP1, y = CAP2, colour = territory)) + geom_point(size = 2)

gg_ordiplot(rda_scores_sib_territory_df_f22_uniq, groups = df_metadata_sib_f22_uniq$territory, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE)

# _ ring ----
rda_sib_ring_f22 = capscale(formula = df_otus_sib_f22 ~ ring_number, data = df_metadata_sib_f22,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_ring_f22) #p = 0.803

# _ age ----
df_metadata_sib_f22 |> group_by(age) |> summarise(n = n()) #adult 19, juvenile 30
rda_sib_age_f22 = capscale(formula = df_otus_sib_f22 ~ age, data = df_metadata_sib_f22,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_age_f22) #p = 0.063
RsquareAdj(rda_sib_age) #adj R^2 = 0.00527

rda_scores_sib_age_df_f22 = rda_sib_age_f22 |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib_f22, by = "study_id")

ggplot(rda_scores_sib_age_df_f22, aes(x = CAP1, y = MDS1, colour = age)) + geom_point()

gg_ordiplot(rda_sib_age, groups = df_metadata_sib$age, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE)

# _ breeding status ----
## all breeders are adults, but non-breeders are both
rda_sib_breeding_f22 = capscale(formula = df_otus_sib_f22 ~ breeding_status, data = df_metadata_sib_f22,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_breeding_f22) #p = 0.278

# _ model selection ----
# on territory and age (since age was approaching significance)
# run null and full models
mod0_sib_f22 = capscale(df_otus_sib_f22 ~ 1, data = df_metadata_sib_f22, distance = "robust.aitchison", na.action = na.exclude)
mod1_sib_f22 = capscale(formula = df_otus_sib_f22 ~ territory + age, data = df_metadata_sib_f22, distance = "robust.aitchison", na.action = na.exclude)

anova(mod1_sib_f22) #p = 0.021
RsquareAdj(mod1_sib_f22) #ajd R^2 = 0.1645092

# model selection on R^2 and p-values
ordiR2step(mod0_sib_f22, scope = formula(mod1_sib_f22), perm.max = 999, na.action = na.exclude) #territory is strongest predictor

# winter 22 vs 23 ----
samples_w22_w23 = df_metadata_sib |> filter(season == "winter2022" | season == "winter2023") |> pull(study_id)
df_otus_w22_w23 = df_otus_sib |> as.data.frame() |> rownames_to_column() |> filter(rowname %in% samples_w22_w23)|> column_to_rownames("rowname")
df_otus_w22_w23 = df_otus_w22_w23[,-(which(colSums(df_otus_w22_w23) == 0))]
df_metadata_w22_w23 = df_metadata_sib |> filter(season == "winter2022"| season == "winter2023")

# _ season ----
rda_sib_winters = capscale(formula = df_otus_w22_w23 ~ season, data = df_metadata_w22_w23,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_winters) #F = 1.0827, p = 0.229
RsquareAdj(rda_sib_winters) #adj R^2 = 0.00343522

# core  ----
# _ setup ----
# set prevalence & detection levels
prevalences = seq(.05, 1, .05)
det = c(0, 0.1, 0.5, 2, 5)/100

# define a common theme for all plots
core_common_theme = theme_classic() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "none")

create_heatmap_plot <- function(pseq_data, prevalences, detections) {
  plot_core(pseq_data,
            plot.type = "heatmap",
            prevalences = prevalences,
            detections = det, 
            min.prevalence = .75)}

customize_plot <- function(plot, title) {
  plot + 
    core_common_theme +  # Apply common theme
    ggtitle(title) +  # Add title
    xlab("Detection threshold") +  # Label x-axis
    ylab("Genus\n") +  # Label y-axis
    scale_fill_viridis() +  # Color scale
    coord_fixed(ratio = 1)  # Adjust aspect ratio
}

# _ genera ----
siberian_rel = microbiome::transform(siberian, "compositional")
siberian_rel_genus = aggregate_taxa(siberian_rel, "Genus")
siberian_rel_genus = subset_taxa(siberian_rel_genus, Genus != "Unknown") #removing unknowns makes it non-compositional
siberian_rel_genus = microbiome::transform(siberian_rel_genus, "compositional")

core_gen_all_seasons <- create_heatmap_plot(siberian_rel_genus, prevalences, det)
core_gen_all_seasons <- customize_plot(core_gen_all_seasons, "All seasons") #(n = 82)
core_gen_all_seasons

# season subsets
core_gen_w22 <- subset_samples(siberian_rel_genus, season == "winter2022") 
core_gen_w22 <- create_heatmap_plot(core_gen_w22, prevalences, det)
core_gen_w22 <- customize_plot(core_gen_w22, "Winter 2022") #(n = 12)
core_gen_w22

core_gen_s22 <- subset_samples(siberian_rel_genus, season == "summer2022") 
core_gen_s22 <- create_heatmap_plot(core_gen_s22, prevalences, det)
core_gen_s22 <- customize_plot(core_gen_s22, "Summer 2022") #(n = 8)
core_gen_s22

core_gen_f22 <- subset_samples(siberian_rel_genus, season == "fall2022") 
core_gen_f22 <- create_heatmap_plot(core_gen_f22, prevalences, det)
core_gen_f22 <- customize_plot(core_gen_f22, "Fall 2022") #(n = 49)
core_gen_f22

core_gen_w23 <- subset_samples(siberian_rel_genus, season == "winter2023") 
core_gen_w23 <- create_heatmap_plot(core_gen_w23, prevalences, det)
core_gen_w23 <- customize_plot(core_gen_w23, "Winter 2023") #(n = 13)
core_gen_w23


# _ species ----
siberian_rel = microbiome::transform(siberian, "compositional")
siberian_rel_sp = aggregate_taxa(siberian_rel, "Species")

core_sp_all_seasons <- create_heatmap_plot(siberian_rel_sp, prevalences, det)
core_sp_all_seasons <- customize_plot(core_sp_all_seasons, "All seasons") #(n = 82)
core_sp_all_seasons

# season subsets
core_sp_w22 <- subset_samples(siberian_rel_sp, season == "winter2022") 
core_sp_w22 <- create_heatmap_plot(core_sp_w22, prevalences, det)
core_sp_w22 <- customize_plot(core_sp_w22, "Winter 2022") #(n = 12)
core_sp_w22

core_sp_s22 <- subset_samples(siberian_rel_sp, season == "summer2022") 
core_sp_s22 <- create_heatmap_plot(core_sp_s22, prevalences, det)
core_sp_s22 <- customize_plot(core_sp_s22, "Summer 2022") #(n = 8)
core_sp_s22

core_sp_f22 <- subset_samples(siberian_rel_sp, season == "fall2022") 
core_sp_f22 <- create_heatmap_plot(core_sp_f22, prevalences, det)
core_sp_f22 <- customize_plot(core_sp_f22, "Fall 2022") #(n = 49)
core_sp_f22

core_sp_w23 <- subset_samples(siberian_rel_sp, season == "winter2023") 
core_sp_w23 <- create_heatmap_plot(core_sp_w23, prevalences, det)
core_sp_w23 <- customize_plot(core_sp_w23, "Winter 2023") #(n = 13)
core_sp_w23

# fix names
core_sp_all_seasons + scale_y_discrete(labels = c("fragi" = "Pseudomonas fragi", "fonticola" = "Serratia fonticola", "antarctica" = "Pseudomonas antarctica", "thermosphacta" = "Brochothrix thermosphacta", "cibarius" = "Psychrobacter cibarius", "coli" = "Escherichia-Shigella coli", "moatsii" = "Mycoplasma moatsii", "inhibens" = "Carnobacterium inhibens", "Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecium" = "Enterococcus faecium", "Bacteria_Firmicutes_Bacilli_Lactobacillales_Streptococcaceae_Lactococcus_lactis" = "Lactococcus lactis"))

# differential abundance ----
# _ f22 vs w23 ----
# subset
f22_w23 = siberian |> 
  subset_samples(season %in% c("fall2022", "winter2023"))

# aldex
v_f22_w23 = df_metadata_sib |> filter(season == "fall2022" | season == "winter2023") |> pull(season) |> as.character()
season_clr = aldex.clr(otu_table(f22_w23), v_f22_w23, mc.samples = 200)
season_ttest = aldex.ttest(season_clr) #takes ~ 1 min
aldex_season_effect = aldex.effect(season_clr, CI = TRUE)
season_aldex_all = data.frame(season_ttest, aldex_season_effect)

# plots
par(mfrow = c(1,2))
aldex.plot(season_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(season_aldex_all, type = "MW", test = "welch", main = "effect plot")

# generate df of differentially abundant ASVs
differential_asvs_season = season_aldex_all |> 
  filter(we.eBH < 0.05) |> 
  rownames_to_column("...1") |> 
  left_join(taxonomy_table, by = "...1") |> # add ASV taxonomical info
  replace_na(list(Genus = "")) |>
  replace_na(list(Species = "sp.")) |>
  mutate(otu_scientific_name = paste(Genus, Species, sep = " ")) |> #concatenate to full ASV name
  mutate(otu_scientific_name = str_trim(otu_scientific_name, side = "left")) |>
  mutate(otu_scientific_name = as.factor(otu_scientific_name)) |> 
  mutate(ci_product = effect.low * effect.high ) |> #column that returns TRUE if CI does not contain 0# |> 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero

# extract results
differential_asvs_season |> 
  mutate(
    summary_text = glue(
      "{otu_scientific_name} was significantly different between conditions ",
      "(effect = {round(effect, 2)}, ",
      "adj. p = {signif(wi.eBH, 3)})."
    )
  ) %>%
  #select(summary_text) |> 
  print()

# plot
differential_asvs_season |> 
  ggplot(aes(x = effect, y = reorder(otu_scientific_name, (effect*-1)))) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
 # ggtitle("Differentially abundant ASVs \nfall 2022 vs winter 2023\n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none",
        #legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"))

# _ w22 vs w23 ----
# subset
w22_w23 = siberian |> 
  subset_samples(season %in% c("winter2022", "winter2023"))

# aldex
v_w22_w23 = df_metadata_sib |> filter(season == "winter2022" | season == "winter2023")|> 
  pull(season) |> as.character()
winter_clr = aldex.clr(otu_table(w22_w23), v_w22_w23, mc.samples = 400)
winter_ttest = aldex.ttest(winter_clr) #takes ~ 1 min
aldex_winter_effect = aldex.effect(winter_clr, CI = TRUE)
winter_aldex_all = data.frame(winter_ttest, aldex_winter_effect)

# plots
par(mfrow = c(1,2))
aldex.plot(winter_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(winter_aldex_all, type = "MW", test = "welch", main = "effect plot")
# there are no differentially abundant ASVs between the two winters
