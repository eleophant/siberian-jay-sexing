###
#
# jays figures
#
###

# packages ----
library(tidyverse)
library(viridis)
library(patchwork)
library(ggmap)
library(ggspatial)
library(sf)
library(ggpubr)
library(cowplot)
library(phyloseq)
library(lme4)
library(vegan)
library(ggordiplots)
library(microbiome)
library(ALDEx2)
theme_set(theme_classic())


# data input ----
df_metadata_sib = read_csv("df_metadata_sib_no_duplicates.csv")
df_otus_sib = read_csv("df_otus_sib_no_duplicates.csv")
taxonomy_table = read_csv("jay_otu_taxonomy_table.csv")
territory_sizes_f22 = read_csv("territories_f22.csv")

# check if any split samples (samples from the same bird on the same day)
df_metadata_sib |> distinct(ring_number, date, .keep_all = TRUE) #there are none

# generate otu matrix
df_otus_sib = df_otus_sib |> column_to_rownames("rowname")

# subset fall22
# subset metadata
df_metadata_sib_f22 = df_metadata_sib |> filter(season == "fall2022")
df_metadata_sib_f22 = df_metadata_sib_f22 |> left_join(territory_sizes_f22, by = join_by(territory))

# subset otus
samples_f22 = df_metadata_sib |> filter(season == "fall2022") |> pull(study_id)
df_otus_sib_f22 = df_otus_sib |> as.data.frame() |> rownames_to_column() |> filter(rowname %in% samples_f22)
df_otus_sib_f22 = df_otus_sib_f22 |> column_to_rownames("rowname")
df_otus_sib_f22 = df_otus_sib_f22[,-(which(colSums(df_otus_sib_f22) == 0))]


# 1) map ----
territories = df_metadata_sib |> 
  select(X, Y) |> 
  rename(longitude = X, latitude = Y) |> 
  distinct(longitude, latitude)

## main map ----
# define the map center
map = get_map(location = map_center, zoom = 10, maptype = "satellite", source = "google")

# create a data frame for site coordinates
sites = data.frame(longitude, latitude)

# calculate reivo vs managed means
north_cluster = sites |>  filter(latitude > 65.74) |>  summarize(
  longitude = mean(longitude),
  latitude = mean(latitude)
)
south_cluster = sites |> filter(latitude <= 65.74) |> summarize(
  longitude = mean(longitude),
  latitude = mean(latitude)
)

# combine clusters into a single data frame with labels
clusters = rbind(
  cbind(north_cluster, label = "protected"),
  cbind(south_cluster, label = "managed")
)

# map
main_map = ggmap(map) +
  # add territory points
  geom_point(data = sites, aes(x = longitude, y = latitude), color = "white", size = 3, alpha = 0.5) +
  geom_text(data = clusters, aes(x = longitude, y = latitude, label = label), 
            nudge_y = 0.04, color = "white", size = 4) +
  
  # move the scale bar to the bottom left
  scale_bar(x = 18.78, y = 65.58, distance_km = 10, text_offset = 0.018, color = "white") +
  
  # add north arrow
  annotation_north_arrow(
    location = "br",  # bottom right corner
    which_north = "true",  # true north
    style = north_arrow_orienteering(line_width = 0, fill = c("white", "white"), text_size = 8, text_col = "white")) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()

main_map


## inset ----
inset_satellite_map = get_map(
  location = c(longitude = 14.6, latitude = 62),  # center of Sweden
  zoom = 5,
  maptype = "satellite",
  source = "google")

inset_map = ggmap(inset_satellite_map) +
  # add a point for the study area
  geom_point(data = data.frame(x = 19.12656, y = 65.78662), 
             aes(x = x, y = y), color = "white", size = 5, alpha = 0.8) +
  theme_void()

inset_map

# 2A) rel ab ----
## normalize number of reads using media sequencing depth
total = median(sample_sums(siberian))
standf = function(x, t=total) round(t * (x / sum(x)))
siberian_norm = transform_sample_counts(siberian, standf)

plot_bar(siberian_norm, fill = "Phylum") +
  geom_bar(aes(colour = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Normalized abundance\n") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12)) +
  ggtitle("A")


# 2B) alpha div obs ----
## season ----
# significantly different seasons
significant_comparisons = list(
  c("winter2022", "fall2022"),
  c("winter2022", "winter2023"))

# calculate sample sizes
df_metadata_sib |> 
  group_by(season) |> 
  summarise(n = n())
#fall2022 49, summer2022 8, winter2022 12, winter2023 13

# set season order
df_metadata_sib = df_metadata_sib |> 
  mutate(season = fct_relevel(season, "winter2022", "summer2022", "fall2022", "winter2023"))

df_metadata_sib |>
  ggplot(aes(x = season, y = Observed, fill = season))+
  geom_boxplot(size = 1) +
  geom_point() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12), axis.text = element_text(size = 12)) +
  xlab("Season") +
  ylab("Observed ASVs") +
  scale_x_discrete(labels=c("winter2022" = "Winter \n2022", "summer2022" = "Summer \n2022", "fall2022" = "Fall \n2022", "winter2023" = "Winter \n2023")) +
  scale_fill_viridis_d(name = "season", labels = c("winter 2022", "summer 2022", "fall 2022", "winter 2023")) +
  
  # add pairwise wilcoxon test
  stat_compare_means(
    method = "wilcox.test",
    comparisons = significant_comparisons,
    label = "p.signif", #use "p.format" for exact p-values
    p.adjust.method = "bonferroni"
  ) +
  ggtitle("B")

## area ----
df_metadata_sib |>
  ggplot(aes(x = area, y = Observed, fill = area))+
  geom_boxplot(size = 1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  scale_fill_viridis_d(name = "area", labels = c("managed", "protected")) +
  xlab("area") +
  ylab("observed ASVs") +
  scale_x_discrete(labels=c("managed" = "managed", "reivo" = "protected"))+
  #add stats
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif", #use "p.format" for exact p-values
  )+
  ggtitle("C")

## age ----
df_metadata_sib |>
  ggplot(aes(x = age, y = Observed, fill = age))+
  geom_boxplot(size = 1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  scale_fill_viridis_d() +
  ylab("observed ASVs") +
  #add stats
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif", #use "p.format" for exact p-values
  ) +
  ggtitle("D")

## breeding status ----
df_metadata_sib |>
  filter(breeding_status != "NA") |> 
  ggplot(aes(x = breeding_status, y = Observed, fill = breeding_status))+
  geom_boxplot(size = 1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  scale_fill_viridis_d(name = "breeding status") +
  xlab("breeding status") +
  ylab("observed ASVs") +
  #scale_x_discrete(labels=c("managed" = "managed", "reivo" = "protected"))+
  #add stats
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif", #use "p.format" for exact p-values
  )+
  ggtitle("E")

# (SI) alpha div shannon ----
## season ----
# significantly different seasons
significant_comparisons = list(
  c("winter2022", "fall2022"),
  c("winter2022", "winter2023"))

df_metadata_sib |>
  ggplot(aes(x = season, y = Shannon, fill = season))+
  geom_boxplot(size = 1) +
  geom_point() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  xlab("season") +
  ylab("Shannon index") +
  scale_x_discrete(labels=c("winter2022" = "winter \n2022", "summer2022" = "summer \n2022", "fall2022" = "fall \n2022", "winter2023" = "winter \n2023")) +
  scale_fill_viridis_d(name = "season", labels = c("winter 2022", "summer 2022", "fall 2022", "winter 2023")) +
  #add pairwise wilcox test
  stat_compare_means(
    method = "wilcox.test",
    comparisons = significant_comparisons,
    label = "p.signif", #use "p.format" for exact p-values
    p.adjust.method = "bonferroni"
  )

## area ----
df_metadata_sib |>
  ggplot(aes(x = area, y = Shannon, fill = area))+
  geom_boxplot(size = 1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  scale_fill_viridis_d(name = "area", labels = c("managed", "protected")) +
  xlab("area") +
  ylab("Shannon index") +
  scale_x_discrete(labels=c("managed" = "managed", "reivo" = "protected"))+
  #add stats
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif", #use "p.format" for exact p-values
  )

## age ----
df_metadata_sib |>
  ggplot(aes(x = age, y = Shannon, fill = age))+
  geom_boxplot(size = 1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  scale_fill_viridis_d() +
  ylab("Shannon index") +
  #add stats
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif", #use "p.format" for exact p-values
  )

## breeding status ----
df_metadata_sib |>
  filter(breeding_status != "NA") |> 
  ggplot(aes(x = breeding_status, y = Shannon, fill = breeding_status))+
  geom_boxplot(size = 1) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  scale_fill_viridis_d(name = "breeding status") +
  xlab("breeding status") +
  ylab("Shannon index") +
  #scale_x_discrete(labels=c("managed" = "managed", "reivo" = "protected"))+
  #add stats
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif", #use "p.format" for exact p-values
  )


# 3A) rda ----
## season ---- 
rda_sib_season = capscale(formula = df_otus_sib ~ season, data = df_metadata_sib,  distance = "robust.aitchison", na.action = na.exclude)
anova(rda_sib_season) #p = 0.001
RsquareAdj(rda_sib_season) #adj R^2 = 0.05109595

rda_scores_sib_season_df = rda_sib_season |> 
  scores(display = "sites") |> as.data.frame() |> rownames_to_column() |> mutate(study_id = rowname) |> 
  left_join(df_metadata_sib, by = "study_id")

# set season order
rda_scores_sib_season_df = rda_scores_sib_season_df |> 
  mutate(season = fct_relevel(season, "winter2022", "summer2022", "fall2022", "winter2023"))

ggplot(rda_scores_sib_season_df, aes(x = CAP1, y = CAP2, colour = season)) + 
  geom_point(size = 2) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d(
    labels = c("winter2022" = "Winter 2022", 
               "summer2022" = "Summer 2022", 
               "fall2022" = "Fall 2022",
               "winter2023" = "Winter 2023")) +
  labs(
    colour = "Season", 
    #shape = "sample type",
    x = "CAP1 (5.74%)", 
    y = "CAP2 (1.68%)"
  ) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("A")


# 3B) diff ab ----
# subset f22 vs w23
f22_w23 = siberian |> 
  subset_samples(season %in% c("fall2022", "winter2023"))

# aldex
v_f22_w23 = df_metadata_sib |> filter(season == "fall2022" | season == "winter2023") |> pull(season) |> as.character()
season_clr = aldex.clr(otu_table(f22_w23), v_f22_w23, mc.samples = 400)
season_ttest = aldex.ttest(season_clr) #takes ~ 1 min
aldex_season_effect = aldex.effect(season_clr, CI = TRUE)
season_aldex_all = data.frame(season_ttest, aldex_season_effect)

# plots
par(mfrow = c(1,2))
aldex.plot(season_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(season_aldex_all, type = "MW", test = "welch", main = "effect plot")

# generate df of differentially abundant ASVs
differential_asvs_season = season_aldex_all |> 
  filter(wi.eBH < 0.05) |> 
  rownames_to_column("...1") |> 
  left_join(taxonomy_table, by = "...1") |> # add ASV taxonomical info
  replace_na(list(Genus = "")) |>
  replace_na(list(Species = "sp.")) |>
  mutate(otu_scientific_name = paste(Genus, Species, sep = " ")) |> #concatenate to full ASV name
  mutate(otu_scientific_name = str_trim(otu_scientific_name, side = "left")) |>
  mutate(otu_scientific_name = as.factor(otu_scientific_name)) |> 
  mutate(ci_product = effect.low * effect.high ) |> #column that returns TRUE if CI does not contain 0# |> 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero

# plot
differential_asvs_season |> 
  ggplot(aes(x = effect, y = reorder(otu_scientific_name, (effect*-1)), colour = Phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  # ggtitle("Differentially abundant ASVs \nfall 2022 vs winter 2023\n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_minimal() +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"),
        aspect.ratio = 1/3) +
  scale_colour_manual(values = c("darkorchid", "darkorange")) +
  ggtitle("B")



# 4) core ----
siberian_rel = microbiome::transform(siberian, "compositional")
siberian_rel_genus = aggregate_taxa(siberian_rel, "Genus")
siberian_rel_genus = subset_taxa(siberian_rel_genus, Genus != "Unknown") #removing unknowns makes it non-compositional
siberian_rel_genus = microbiome::transform(siberian_rel_genus, "compositional")

# set prevalence & detection levels
prevalences = seq(.05, 1, .05)
det = c(0, 0.1, 0.5, 2, 5)/100

# define a common theme for all plots
core_common_theme = theme_classic() +
  theme(axis.text.y = element_text(face = "italic"))

create_heatmap_plot = function(pseq_data, prevalences, detections) {
  plot_core(pseq_data,
            plot.type = "heatmap",
            prevalences = prevalences,
            detections = det, 
            min.prevalence = .75)}

customize_plot = function(plot, title) {
  plot + 
    core_common_theme + 
    ggtitle(title) +
    xlab("Detection threshold") +
    ylab("Genus\n") +
    scale_fill_viridis_c(option = "inferno") +
    theme(title = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    coord_fixed(ratio = 1)  #adjust aspect ratio
}

## all seasons ----
core_gen_all_seasons = create_heatmap_plot(siberian_rel_genus, prevalences, det)
core_gen_all_seasons = customize_plot(core_gen_all_seasons, "All seasons") #(n = 82)
core_gen_all_seasons

## by season ----
core_gen_w22 = subset_samples(siberian_rel_genus, season == "winter2022") 
core_gen_w22 = create_heatmap_plot(core_gen_w22, prevalences, det)
core_gen_w22 = customize_plot(core_gen_w22, "Winter 2022") #(n = 12)
core_gen_w22

core_gen_s22 = subset_samples(siberian_rel_genus, season == "summer2022") 
core_gen_s22 = create_heatmap_plot(core_gen_s22, prevalences, det)
core_gen_s22 = customize_plot(core_gen_s22, "Summer 2022") #(n = 8)
core_gen_s22

core_gen_f22 = subset_samples(siberian_rel_genus, season == "fall2022") 
core_gen_f22 = create_heatmap_plot(core_gen_f22, prevalences, det)
core_gen_f22 = customize_plot(core_gen_f22, "Fall 2022") #(n = 49)
core_gen_f22

core_gen_w23 = subset_samples(siberian_rel_genus, season == "fall2022") 
core_gen_w23 = create_heatmap_plot(core_gen_w23, prevalences, det)
core_gen_w23 = customize_plot(core_gen_w23, "Winter 2023") #(n = 13)
core_gen_w23

core_all = (core_gen_w22 | core_gen_s22) /
  (core_gen_f22 | core_gen_w23)
core_all
