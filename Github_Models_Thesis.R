library(sf)
library(dplyr)
library(tidyr)

names(trails_disturb)
head(trails_disturb)

#create camera points spatial object
camera_points <- master_presence %>%
  distinct(site_id, .keep_all = TRUE) %>% #1RET023 uses just one of these since both deployements with this same ID was have same coords, diff when creating a matrix where both capture diff animals.
  select(site_id, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(trails_disturb))  # match trail CRS

#define buffer sizes
buffer_sizes <- c(100, 250, 500, 1000)

#calculate trail length within each buffer
trail_buffer_results <- list()

#create buffers around cameras
for(buf in buffer_sizes) {camera_buffers <- camera_points %>%
  st_buffer(dist = buf) %>% select(site_id)

#intersect buffers with trails
trails_clipped <- st_intersection(trails_disturb, camera_buffers)

#calculate length of trail within each buffer
#split by high/low use
trail_lengths <- trails_clipped %>%
  mutate(trail_length_m = as.numeric(st_length(.))) %>%
  group_by(site_id, Disturb) %>%
  summarise(trail_length_m = sum(trail_length_m, na.rm = TRUE),
            .groups = "drop") %>%
  st_drop_geometry() %>%
  pivot_wider(names_from  = Disturb,
              values_from = trail_length_m,
              values_fill = 0) %>%
  rename(high_use_trail_m = High,
         low_use_trail_m  = Low) %>%
  mutate(total_trail_m = high_use_trail_m + low_use_trail_m,
         buffer_m      = buf)

#add cameras with zero trail length
all_cameras <- camera_points %>%
  st_drop_geometry() %>%
  select(site_id)

trail_lengths_complete <- all_cameras %>%
  left_join(trail_lengths, by = "site_id") %>%
  mutate(
    buffer_m         = buf,
    high_use_trail_m = replace_na(high_use_trail_m, 0),
    low_use_trail_m  = replace_na(low_use_trail_m, 0),
    total_trail_m    = replace_na(total_trail_m, 0))

trail_buffer_results[[as.character(buf)]] <- trail_lengths_complete}

#combine into one table
trail_buffer_table <- bind_rows(trail_buffer_results) %>%
  select(site_id, buffer_m, high_use_trail_m, 
         low_use_trail_m, total_trail_m) %>%
  arrange(site_id, buffer_m)

#check
head(trail_buffer_table, 16)
nrow(trail_buffer_table)  #should be 104 sites x 4 buffers = 416 rows

#checking distribution of lengths, if there are too many 0s at 100m then we can skip 100m
trail_buffer_table %>%
  group_by(buffer_m) %>%
  summarise(sites_with_trail = sum(total_trail_m > 0),
            sites_with_high  = sum(high_use_trail_m > 0),
            median_total     = median(total_trail_m))
#buffer_m sites_with_trail sites_with_high median_total
#<dbl>            <int>           <int>        <dbl>
#1      100               18               9           0 
#2      250               33              19           0 
#3      500               58              37         275.
#4     1000               82              52        2302.

#checking which buffer size overlaps, with cameras near by. 
camera_buffers_500 <- st_buffer(camera_points, dist = 500)
overlaps_500 <- st_intersects(camera_buffers_500, camera_buffers_500)
sum(lengths(overlaps_500) > 1)   # cameras overlapping at 500m

camera_buffers_1000 <- st_buffer(camera_points, dist = 1000)
overlaps_1000 <- st_intersects(camera_buffers_1000, camera_buffers_1000)
sum(lengths(overlaps_1000) > 1)  # cameras overlapping at 1000m

#checking distribution
trail_500 <- trail_buffer_table %>% filter(buffer_m == 500)
summary(trail_500$total_trail_m)
hist(trail_500$total_trail_m, breaks = 20)

trail_500 %>%
  left_join(trail_covariates %>% select(site_id, distance_to_trail_m),
            by = "site_id") %>%
  summarise(r = cor(total_trail_m, log1p(distance_to_trail_m)))

hist(log1p(trail_500$total_trail_m), breaks = 20)

##are we log transformaing or no? I think no because values get all funky, but will scale.

trail_500 <- trail_500 %>%
  mutate(trail_500_sc = as.numeric(scale(total_trail_m)))

site_covs <- site_covs %>%
  left_join(trail_500 %>% select(site_id, total_trail_m, trail_500_sc),
            by = "site_id")

summary(site_covs$trail_500_sc)


###### REBUILDING 02 AND 03 with new disturbance variable ############################################

##BROWN BEAR SPATIAL
site_covs_umf <- site_covs %>%
  select(habitat_simple, sheep_rate_sc,
         log_distance_trail_sc, trail_500_sc, elevation_sc)

umf_bear <- unmarkedFramePCount(
  y        = as.matrix(bear_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_bear))

#just checking if the p values are more significant distance vs buffer
nm_dist_bear <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                         log_distance_trail_sc, data = umf_bear, K = 100)
nm_trail500_bear <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                             trail_500_sc, data = umf_bear, K = 100)
summary(nm_trail500_bear)


nm0_bear <- pcount(~ log_effort ~ 1, data = umf_bear, K = 100)
nm_envr_bear <- pcount(~ log_effort ~ habitat_simple + elevation_sc, data = umf_bear, K = 100)
nm_pastoral_bear <- pcount(~ log_effort ~ sheep_rate_sc, data = umf_bear, K = 100)
nm_recreation_bear <- pcount(~ log_effort ~ trail_500_sc, data = umf_bear, K = 100)
nm_envr_pastoral_bear <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                  sheep_rate_sc, data = umf_bear, K = 100)
nm_envr_recreation_bear <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                    trail_500_sc, data = umf_bear, K = 100)
nm_recreation_pastoral_bear <- pcount(~ log_effort ~ trail_500_sc +
                                        sheep_rate_sc, data = umf_bear, K = 100)
nm_full_bear <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                         trail_500_sc + sheep_rate_sc, data = umf_bear, K = 100)

nm_model_list_bear <- fitList(
  "null"                = nm0_bear,
  "envr"                = nm_envr_bear,
  "pastoral"            = nm_pastoral_bear,
  "recreation"          = nm_recreation_bear,
  "envr+pastoral"       = nm_envr_pastoral_bear,
  "envr+recreation"     = nm_envr_recreation_bear,
  "recreation+pastoral" = nm_recreation_pastoral_bear,
  "full"                = nm_full_bear)

nm_model_table_bear <- modSel(nm_model_list_bear)
print(nm_model_table_bear)
#                    nPars     AIC delta   AICwt cumltvWt
#full                    9 1186.23  0.00 6.6e-01     0.66 only Δ1.34 difference between envr+rec, pick simpler model. 
#envr+recreation         8 1187.57  1.34 3.4e-01     1.00 #best
#envr+pastoral           8 1198.45 12.22 1.5e-03     1.00
#envr                    7 1198.62 12.39 1.3e-03     1.00
#recreation              4 1206.29 20.06 2.9e-05     1.00
#recreation+pastoral     5 1208.22 21.98 1.1e-05     1.00
#null                    3 1214.15 27.91 5.7e-07     1.00
#pastoral                4 1215.97 29.74 2.3e-07     1.00
summary(nm_full_bear)
#Abundance (log-scale):
#Estimate    SE        z  P(>|z|)
#(Intercept)               1.217198 0.224  5.43360 5.52e-08
#habitat_simpleRocks      -1.123664 0.358 -3.14008 1.69e-03
#habitat_simpleDwarf_pine -0.000803 0.253 -0.00317 9.97e-01
#habitat_simpleForest     -0.855450 0.263 -3.24873 1.16e-03
#elevation_sc              0.136471 0.116  1.17333 2.41e-01
#trail_500_sc             -0.408437 0.117 -3.50215 4.62e-04 #0.000462
#sheep_rate_sc            -0.191908 0.116 -1.64930 9.91e-02 #0.0991 - sheep not adding anything sig, like in AIC, so best is envr+rec.
summary(nm_envr_recreation_bear)
# Call:
#   pcount(formula = ~log_effort ~ habitat_simple + elevation_sc + 
#            trail_500_sc, data = umf_bear, K = 100)
# 
# Abundance (log-scale):
#   Estimate    SE      z  P(>|z|)
# (Intercept)                 1.069 0.214  5.004 5.62e-07
# habitat_simpleRocks        -0.949 0.350 -2.712 6.68e-03
# habitat_simpleDwarf_pine    0.146 0.246  0.594 5.52e-01
# habitat_simpleForest       -0.667 0.249 -2.677 7.43e-03
# elevation_sc                0.150 0.120  1.244 2.13e-01
# trail_500_sc               -0.393 0.117 -3.361 7.76e-04
# 
# Detection (logit-scale):
#   Estimate   SE     z  P(>|z|)
# (Intercept)    -7.41 1.95 -3.81 0.000139
# log_effort      2.55 1.00  2.55 0.010805
# 
# AIC: 1187.574 
# Number of sites: 104


##BROWN BEAR TEMPORAL
library(lubridate)
library(dplyr)
library(tidyr)
library(lubridate)
library(GLMMadaptive)
source("https://raw.githubusercontent.com/MarcusRowcliffe/make_chm_data/refs/heads/main/make_chm_data.R")

#adding trail density (500m) to deployment covariates
#trail_500 is keyed on site_id matches locationName here?
deployment_covariates_chm <- deployment_covariates_chm %>%
  left_join(trail_500 %>% select(site_id, total_trail_m),
            by = c("locationName" = "site_id"))

#checking NAs
summary(deployment_covariates_chm$total_trail_m)
sum(is.na(deployment_covariates_chm$total_trail_m))

#removing 8 minute deployement
deployment_covariates_chm <- deployment_covariates_chm %>%
  filter(locationName != "2RET071")
sum(is.na(deployment_covariates_chm$total_trail_m))   # should be 0
length(unique(deployment_covariates_chm$locationName)) # should be 104

##habitat-controlled, trail_500_sc replacing human_rate_sc
chm_data_bear <- make_chm_data(
  deployments  = deployment_covariates_chm,
  observations = observations_chm %>%
    filter(scientificName == "Ursus arctos"), nBins = 24,
  covs = c("sheep_rate_100tn", "total_trail_m",
           "habitat_simple"),collapse = TRUE) %>%
  mutate(
    sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
    trail_500_sc  = as.numeric(scale(total_trail_m)),
    habitat_simple = case_when(
      locationName == "1RET022" ~ "Open_pasture",
      locationName == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple)) %>% 
      factor(levels = c("Open_pasture", "Rocks","Dwarf_pine", "Forest")))

#checks
length(unique(chm_data_bear$locationName))   #104
table(chm_data_bear$habitat_simple)
summary(chm_data_bear$trail_500_sc)

#Step 1 - base models
chm_null_bear <- mixed_model(
  fixed  = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_bear)

chm_unimodal_bear <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_bear)

chm_bimodal_bear <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_bear)

#Step 2 - habitat added to bimodal
chm_bimodal_hab_bear <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_bear)

#Step 3 - disturbance controlling for habitat
chm_sheep_hab_bear <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_bear)

#NEW - trail density replaces human rate
chm_trail_hab_bear <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_bear)

chm_trail_sheep_hab_bear <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_bear)

AIC(chm_null_bear, chm_unimodal_bear, chm_bimodal_bear,
    chm_bimodal_hab_bear, chm_sheep_hab_bear, chm_trail_hab_bear, chm_trail_sheep_hab_bear)
#df      AIC
#chm_null_bear             2 1290.491
#chm_unimodal_bear         4 1267.049
#chm_bimodal_bear          6 1239.708
#chm_bimodal_hab_bear      9 1235.142 #1235.14 → 1221.39 = Δ13.75
#chm_sheep_hab_bear       14 1221.387 #best
#chm_trail_hab_bear       14 1235.600
#chm_trail_sheep_hab_bear 19 1220.458 #best but not by much compared to simple sheep. 
summary(chm_sheep_hab_bear)
#                                 Estimate Std.Err  z-value  p-value
#(Intercept)                        -7.6828  0.3264 -23.5386  < 1e-04
#cos(timeRadian)                    -0.3858  0.1555  -2.4819 0.013068
#sheep_rate_sc                      -0.8591  0.3540  -2.4264 0.015248 (main effect)
#sin(timeRadian)                    -0.2348  0.1048  -2.2411 0.025020
#cos(2 * timeRadian)                -0.6648  0.1266  -5.2491  < 1e-04
#sin(2 * timeRadian)                 0.1181  0.1260   0.9377 0.348408
#habitat_simpleRocks                -0.8108  0.5338  -1.5189 0.128780
#habitat_simpleDwarf_pine            0.6397  0.4422   1.4465 0.148031
#habitat_simpleForest               -0.6595  0.3930  -1.6782 0.093301
#cos(timeRadian):sheep_rate_sc       0.6706  0.4029   1.6645 0.096013
#sheep_rate_sc:sin(timeRadian)       0.0212  0.2328   0.0910 0.927523
#sheep_rate_sc:cos(2 * timeRadian)  -0.3686  0.2788  -1.3219 0.186200
#sheep_rate_sc:sin(2 * timeRadian)  -0.8144  0.3229  -2.5221 0.011667 (a significant interaction?)
summary(chm_trail_sheep_hab_bear)
# Estimate Std.Err  z-value   p-value
# (Intercept)                        -7.6146  0.3257 -23.3759   < 1e-04
# cos(timeRadian)                    -0.3715  0.1698  -2.1885 0.0286298
# trail_500_sc                       -0.4732  0.1766  -2.6793 0.0073777
# sin(timeRadian)                    -0.2437  0.1079  -2.2584 0.0239180
# cos(2 * timeRadian)                -0.7605  0.1384  -5.4959   < 1e-04
# sin(2 * timeRadian)                 0.1223  0.1335   0.9160 0.3596434
# sheep_rate_sc                      -0.9025  0.3561  -2.5344 0.0112639
# habitat_simpleRocks                -0.9878  0.5314  -1.8589 0.0630447
# habitat_simpleDwarf_pine            0.4607  0.4380   1.0519 0.2928660
# habitat_simpleForest               -0.8252  0.3943  -2.0925 0.0363976
# cos(timeRadian):trail_500_sc        0.0341  0.1776   0.1921 0.8476794
# trail_500_sc:sin(timeRadian)       -0.0549  0.1102  -0.4986 0.6180663
# trail_500_sc:cos(2 * timeRadian)   -0.3049  0.1511  -2.0176 0.0436334
# trail_500_sc:sin(2 * timeRadian)    0.0012  0.1358   0.0086 0.9931351
# cos(timeRadian):sheep_rate_sc       0.7014  0.4009   1.7496 0.0801838
# sin(timeRadian):sheep_rate_sc       0.0180  0.2345   0.0769 0.9386714
# cos(2 * timeRadian):sheep_rate_sc  -0.3812  0.2805  -1.3590 0.1741569
# sin(2 * timeRadian):sheep_rate_sc  -0.8315  0.3228  -2.5760 0.0099963

#fit_chm comparason
fitchm_trail_habint_bear <- fit_chm(
  cbind(success, failure) ~ trail_500_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_bear)

fitchm_sheep_habint_bear <- fit_chm(
  cbind(success, failure) ~ sheep_rate_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_bear)

AIC(chm_bimodal_hab_bear,
    chm_trail_hab_bear,
    chm_sheep_hab_bear,
    fitchm_trail_habint_bear,
    fitchm_sheep_habint_bear)
summary(fitchm_sheep_habint_bear)

##RED DEER SPATIAL
umf_rd <- unmarkedFramePCount(
  y = as.matrix(red_deer_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_rd))

nm0_rd <- pcount(~ log_effort ~ 1, data = umf_rd, K = 100)
nm_envr_rd <- pcount(~ log_effort ~ habitat_simple + elevation_sc,
                     data = umf_rd, K = 100)
nm_pastoral_rd <- pcount(~ log_effort ~ sheep_rate_sc,
                         data = umf_rd, K = 100)
nm_recreation_rd <- pcount(~ log_effort ~ trail_500_sc,
                           data = umf_rd, K = 100)
nm_envr_pastoral_rd <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                sheep_rate_sc, data = umf_rd, K = 100)
nm_envr_recreation_rd <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                  trail_500_sc, data = umf_rd, K = 100)
nm_recreation_pastoral_rd <- pcount(~ log_effort ~ trail_500_sc +
                                      sheep_rate_sc, data = umf_rd, K = 100)
nm_full_rd <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                       trail_500_sc + sheep_rate_sc, data = umf_rd, K = 100)
nm_model_list_rd <- fitList(
  "null"                = nm0_rd,
  "envr"                = nm_envr_rd,
  "pastoral"            = nm_pastoral_rd,
  "recreation"          = nm_recreation_rd,
  "envr+pastoral"       = nm_envr_pastoral_rd,
  "envr+recreation"     = nm_envr_recreation_rd,
  "recreation+pastoral" = nm_recreation_pastoral_rd,
  "full"                = nm_full_rd)

nm_model_table_rd <- modSel(nm_model_list_rd)
print(nm_model_table_rd)
#                    nPars     AIC delta   AICwt cumltvWt
#envr+recreation         8 1649.58  0.00 6.0e-01     0.60 #best
#full                    9 1650.62  1.04 3.6e-01     0.96 #Arnold thing here too, similar sheep not being sig like in bear sp.
#envr                    7 1656.17  6.59 2.2e-02     0.99
#envr+pastoral           8 1657.20  7.62 1.3e-02     1.00
#recreation              4 1682.90 33.32 3.5e-08     1.00
#recreation+pastoral     5 1684.46 34.87 1.6e-08     1.00
#null                    3 1687.77 38.19 3.1e-09     1.00
#pastoral                4 1689.57 39.99 1.3e-09     1.00
summary(nm_full_rd)
summary(nm_envr_recreation_rd)
# Call:
#   pcount(formula = ~log_effort ~ habitat_simple + elevation_sc + 
#            trail_500_sc, data = umf_rd, K = 100)
# 
# Abundance (log-scale):
#   Estimate    SE     z  P(>|z|)
# (Intercept)                 0.471 0.197  2.39 0.017068
# habitat_simpleRocks        -2.100 0.609 -3.45 0.000561
# habitat_simpleDwarf_pine    0.365 0.256  1.42 0.154283
# habitat_simpleForest        0.468 0.235  1.99 0.046238
# elevation_sc                0.258 0.106  2.45 0.014444
# trail_500_sc               -0.267 0.096 -2.79 0.005345
# 
# Detection (logit-scale):
#   Estimate    SE     z  P(>|z|)
# (Intercept)    -3.33 0.334 -9.95 2.58e-23
# log_effort      0.81 0.169  4.78 1.75e-06
# 
# AIC: 1649.583 
# Number of sites: 104

##RED DEER TEMPORAL
chm_data_red_deer <- make_chm_data(deployments  = deployment_covariates_chm,
                                   observations = observations_chm %>%
                                     filter(scientificName == "Cervus elaphus"),
                                   nBins = 24,
                                   covs  = c("sheep_rate_100tn", "total_trail_m",
                                             "habitat_simple"),
                                   collapse = TRUE) %>% mutate(
                                     sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
                                     trail_500_sc  = as.numeric(scale(total_trail_m)),
                                     habitat_simple = case_when(
                                       locationName == "1RET022" ~ "Open_pasture",
                                       locationName == "1RET035" ~ "Dwarf_pine",
                                       TRUE ~ as.character(habitat_simple)) %>% 
                                       factor(levels = c("Open_pasture", "Rocks", "Dwarf_pine", "Forest")))

length(unique(chm_data_red_deer$locationName))   #104

#base mods
chm_null_rd <- mixed_model(
  fixed  = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_red_deer)

chm_unimodal_rd <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_red_deer)

chm_bimodal_rd <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_red_deer)

#habitat control
chm_bimodal_hab_rd <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_red_deer)

#disturbance, habitat controlled
chm_sheep_hab_rd <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_red_deer)

chm_trail_hab_rd <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_red_deer)

chm_trail_sheep_hab_rd <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_red_deer)

AIC(chm_null_rd, chm_unimodal_rd, chm_bimodal_rd,
    chm_bimodal_hab_rd, chm_sheep_hab_rd, chm_trail_hab_rd, chm_trail_sheep_hab_rd)
# df      AIC
# chm_null_rd             2 1581.217
# chm_unimodal_rd         4 1564.657
# chm_bimodal_rd          6 1537.940
# chm_bimodal_hab_rd      9 1528.534 #1528.53 → 1517.72 = Δ10.81
# chm_sheep_hab_rd       14 1517.720 #best
# chm_trail_hab_rd       14 1534.285
# chm_trail_sheep_hab_rd 19 1523.688
summary(chm_sheep_hab_rd)
#                                  Estimate Std.Err  z-value  p-value
# (Intercept)                        -8.1709  0.4773 -17.1190  < 1e-04
# cos(timeRadian)                     0.5210  0.1106   4.7129  < 1e-04
# sheep_rate_sc                      -0.2165  0.2646  -0.8182 0.413249
# sin(timeRadian)                     0.0289  0.0788   0.3674 0.713311
# cos(2 * timeRadian)                -0.4962  0.0903  -5.4945  < 1e-04
# sin(2 * timeRadian)                 0.1256  0.0864   1.4528 0.146273
# habitat_simpleRocks                -2.2790  0.9631  -2.3663 0.017968
# habitat_simpleDwarf_pine            0.8934  0.6463   1.3824 0.166846
# habitat_simpleForest                0.5898  0.5519   1.0687 0.285214
# cos(timeRadian):sheep_rate_sc       0.5730  0.2461   2.3280 0.019911
# sheep_rate_sc:sin(timeRadian)       0.0720  0.1225   0.5878 0.556686
# sheep_rate_sc:cos(2 * timeRadian)  -0.1070  0.1378  -0.7759 0.437824
# sheep_rate_sc:sin(2 * timeRadian)  -0.2436  0.1188  -2.0503 0.040338

#fit_chm comparason
fitchm_trail_habint_rd <- fit_chm(
  cbind(success, failure) ~ trail_500_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_red_deer)
# Error in mixed_fit(y, X, Z, X_zi, Z_zi, id, offset, offset_zi, family,  : 
# A large coefficient value has been detected during the optimization. Please
# re-scale you covariates and/or try setting the control argument 'iter_EM = 0'.
# Alternatively, this may due to a divergence of the optimization algorithm,
# indicating that an overly complex model is fitted to the data. For example,
# this could be caused when including random-effects terms (e.g., in the
# zero-inflated part) that you do not need. Otherwise, adjust the
# 'max_coef_value' control argument.
#rocky hab only has 3 rd detecs. with extra parameters in fit_chm model fails. 
#it works again when we drop rocky habitat.
# habitat_simple n_bins successes bins_with_detection
# 1 Open_pasture      576        65                  47
# 2 Rocks             336         3                   3
# 3 Dwarf_pine        408        73                  58
# 4 Forest           1176       164                 122
#Drop rocks, refit
#Drop rocks, refit
fitchm_trail_habint_rd_test <- fit_chm(
  cbind(success, failure) ~ trail_500_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_red_deer %>% filter(habitat_simple != "Rocks") %>% droplevels())

fitchm_sheep_habint_rd_test <- fit_chm(
  cbind(success, failure) ~ sheep_rate_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_red_deer %>%
    filter(habitat_simple != "Rocks") %>%
    droplevels())

AIC(chm_bimodal_hab_rd,
    chm_trail_hab_rd,
    chm_sheep_hab_rd,
    fitchm_trail_habint_rd,
    fitchm_sheep_habint_rd)


##RED FOX SPATIAL
umf_red_fox <- unmarkedFramePCount(
  y = as.matrix(red_fox_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_red_fox))

nm0_red_fox <- pcount(~ log_effort ~ 1, data = umf_red_fox, K = 100)
nm_envr_red_fox <- pcount(~ log_effort ~ habitat_simple + elevation_sc,
                          data = umf_red_fox, K = 100)
nm_pastoral_red_fox <- pcount(~ log_effort ~ sheep_rate_sc,
                              data = umf_red_fox, K = 100)
nm_recreation_red_fox <- pcount(~ log_effort ~ trail_500_sc,
                                data = umf_red_fox, K = 100)
nm_envr_pastoral_red_fox <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                     sheep_rate_sc, data = umf_red_fox, K = 100)
nm_envr_recreation_red_fox <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                       trail_500_sc, data = umf_red_fox, K = 100)
nm_recreation_pastoral_red_fox <- pcount(~ log_effort ~ trail_500_sc +
                                           sheep_rate_sc, data = umf_red_fox, K = 100)
nm_full_red_fox <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                            trail_500_sc + sheep_rate_sc, data = umf_red_fox, K = 100)

nm_model_list_red_fox <- fitList(
  "null"                = nm0_red_fox,
  "envr"                = nm_envr_red_fox,
  "pastoral"            = nm_pastoral_red_fox,
  "recreation"          = nm_recreation_red_fox,
  "envr+pastoral"       = nm_envr_pastoral_red_fox,
  "envr+recreation"     = nm_envr_recreation_red_fox,
  "recreation+pastoral" = nm_recreation_pastoral_red_fox,
  "full"                = nm_full_red_fox)

nm_model_table_red_fox <- modSel(nm_model_list_red_fox)
print(nm_model_table_red_fox)
#                    nPars     AIC delta   AICwt cumltvWt
#envr+pastoral           8 1067.98  0.00 0.38015     0.38
#envr                    7 1068.50  0.52 0.29276     0.67
#full                    9 1069.57  1.59 0.17137     0.84
#envr+recreation         8 1070.05  2.07 0.13492     0.98
#pastoral                4 1074.53  6.55 0.01439     0.99
#recreation+pastoral     5 1076.53  8.55 0.00529     1.00
#null                    3 1080.28 12.30 0.00081     1.00
#recreation              4 1082.23 14.24 0.00031     1.00
summary(nm_envr_pastoral_red_fox)
# Call:
#   pcount(formula = ~log_effort ~ habitat_simple + elevation_sc + 
#            sheep_rate_sc, data = umf_red_fox, K = 100)
# 
# Abundance (log-scale):
#   Estimate     SE        z P(>|z|)
# (Intercept)               0.45380 0.2280  1.99038  0.0465
# habitat_simpleRocks      -0.90389 0.4131 -2.18793  0.0287
# habitat_simpleDwarf_pine -0.00216 0.2975 -0.00727  0.9942
# habitat_simpleForest     -0.67778 0.3120 -2.17231  0.0298
# elevation_sc              0.15181 0.1398  1.08597  0.2775
# sheep_rate_sc             0.12351 0.0723  1.70851  0.0875
# 
# Detection (logit-scale):
#   Estimate    SE     z  P(>|z|)
# (Intercept)    -4.05 0.660 -6.14 8.23e-10
# log_effort      1.20 0.339  3.56 3.78e-04
# 
# AIC: 1067.981 
# Number of sites: 104
summary(nm_envr_red_fox)
# Call:
#   pcount(formula = ~log_effort ~ habitat_simple + elevation_sc, 
#          data = umf_red_fox, K = 100)
# 
# Abundance (log-scale):
#   Estimate    SE      z P(>|z|)
# (Intercept)                 0.581 0.204  2.845 0.00444
# habitat_simpleRocks        -1.054 0.397 -2.653 0.00797
# habitat_simpleDwarf_pine   -0.141 0.278 -0.505 0.61346
# habitat_simpleForest       -0.847 0.287 -2.952 0.00316
# elevation_sc                0.136 0.135  1.006 0.31439
# 
# Detection (logit-scale):
#   Estimate    SE     z  P(>|z|)
# (Intercept)    -4.05 0.660 -6.14 8.43e-10
# log_effort      1.21 0.339  3.56 3.73e-04
# 
# AIC: 1068.504 
# Number of sites: 104

##RED FOX TEMPORAL
chm_data_fox <- make_chm_data(
  deployments  = deployment_covariates_chm,
  observations = observations_chm %>%
    filter(scientificName == "Vulpes vulpes"),
  nBins = 24,
  covs = c("sheep_rate_100tn", "total_trail_m",
           "habitat_simple"),
  collapse = TRUE) %>% mutate(
    sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
    trail_500_sc  = as.numeric(scale(total_trail_m)),
    habitat_simple = case_when(
      locationName == "1RET022" ~ "Open_pasture",
      locationName == "1RET035" ~ "Dwarf_pine",
      TRUE ~ as.character(habitat_simple)) %>% 
      factor(levels = c("Open_pasture", "Rocks",
                        "Dwarf_pine", "Forest")))

length(unique(chm_data_fox$locationName))   #104

#base models
chm_null_fox <- mixed_model(
  fixed  = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_fox)

chm_unimodal_fox <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_fox)

chm_bimodal_fox <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_fox)

#habitat control
chm_bimodal_hab_fox <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_fox)

#disturbance, habitat controlled
chm_sheep_hab_fox <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_fox)

chm_trail_hab_fox <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_fox)

chm_trail_sheep_hab_fox <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_fox)

AIC(chm_null_fox, chm_unimodal_fox, chm_bimodal_fox,
    chm_bimodal_hab_fox, chm_sheep_hab_fox, chm_trail_hab_fox, chm_trail_sheep_hab_fox)
#                    df      AIC
#chm_null_fox         2 1229.568
#chm_unimodal_fox     4 1186.452
#chm_bimodal_fox      6 1187.898
#chm_bimodal_hab_fox  9 1181.018 #1181.02 → 1178.93 = Δ2.09 #not that much
#chm_sheep_hab_fox   14 1178.928 #best
#chm_trail_hab_fox   14 1186.849
#chm_trail_sheep_hab_fox 19 1184.580
summary(chm_sheep_hab_fox)
#                                  Estimate Std.Err  z-value   p-value
#(Intercept)                        -7.6512  0.4327 -17.6822   < 1e-04
# cos(timeRadian)                     0.5819  0.1172   4.9629   < 1e-04
# sheep_rate_sc                       0.1573  0.1988   0.7910 0.4289486
# sin(timeRadian)                    -0.3519  0.1142  -3.0807 0.0020653
# cos(2 * timeRadian)                -0.0087  0.1094  -0.0791 0.9369217
# sin(2 * timeRadian)                -0.1996  0.1094  -1.8248 0.0680383
# habitat_simpleRocks                -1.1835  0.7306  -1.6198 0.1052769
# habitat_simpleDwarf_pine            0.2249  0.6118   0.3675 0.7132365
# habitat_simpleForest               -1.2253  0.5445  -2.2503 0.0244299
# cos(timeRadian):sheep_rate_sc       0.0814  0.1129   0.7208 0.4710453
# sheep_rate_sc:sin(timeRadian)       0.2559  0.1099   2.3278 0.0199232
# sheep_rate_sc:cos(2 * timeRadian)   0.0083  0.0951   0.0874 0.9303621
# sheep_rate_sc:sin(2 * timeRadian)   0.0687  0.0958   0.7180 0.4727648

#fit_chm comparason
fitchm_sheep_habint_fox <- fit_chm(
  cbind(success, failure) ~ sheep_rate_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_fox)

fitchm_trail_habint_fox <- fit_chm(
  cbind(success, failure) ~ trail_500_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_fox)

AIC(chm_bimodal_hab_fox, chm_sheep_hab_fox, chm_trail_hab_fox,
    fitchm_sheep_habint_fox, fitchm_trail_habint_fox)
# df      AIC
# chm_bimodal_hab_fox      9 1181.018
# chm_sheep_hab_fox       14 1178.928
# chm_trail_hab_fox       14 1186.849
# fitchm_sheep_habint_fox 26 1191.741
# fitchm_trail_habint_fox 26 1198.066

##CHAMOIS SPATIAL
umf_chamois_count <- unmarkedFramePCount(
  y  = as.matrix(chamois_count_matrix),
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort))

nm0 <- pcount(~ log_effort ~ 1, data = umf_chamois_count, K = 100)
nm_envr <- pcount(~ log_effort ~ habitat_simple + elevation_sc,
                  data = umf_chamois_count, K = 100)
nm_pastoral <- pcount(~ log_effort ~ sheep_rate_sc,
                      data = umf_chamois_count, K = 100)
nm_recreation <- pcount(~ log_effort ~ trail_500_sc,
                        data = umf_chamois_count, K = 100)
nm_envr_pastoral <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                             sheep_rate_sc, data = umf_chamois_count, K = 100)
nm_envr_recreation <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                               trail_500_sc, data = umf_chamois_count, K = 100)
nm_recreation_pastoral <- pcount(~ log_effort ~ trail_500_sc +
                                   sheep_rate_sc, data = umf_chamois_count, K = 100)
nm_full <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                    trail_500_sc + sheep_rate_sc, data = umf_chamois_count, K = 100)

nm_model_list_chamois <- fitList(
  "null"                = nm0,
  "envr"                = nm_envr,
  "pastoral"            = nm_pastoral,
  "recreation"          = nm_recreation,
  "envr+pastoral"       = nm_envr_pastoral,
  "envr+recreation"     = nm_envr_recreation,
  "recreation+pastoral" = nm_recreation_pastoral,
  "full"                = nm_full)

nm_model_table_chamois <- modSel(nm_model_list_chamois)
print(nm_model_table_chamois)
#                    nPars     AIC  delta   AICwt cumltvWt
#envr+recreation         8 1922.96   0.00 6.6e-01     0.66
#full                    9 1924.26   1.30 3.4e-01     1.00
#envr                    7 1941.46  18.50 6.3e-05     1.00
#envr+pastoral           8 1943.08  20.12 2.8e-05     1.00
#recreation+pastoral     5 2283.59 360.63 3.2e-79     1.00
#recreation              4 2286.68 363.72 6.9e-80     1.00
#pastoral                4 2315.29 392.33 4.2e-86     1.00
#null                    3 2319.13 396.17 6.2e-87     1.00
summary(nm_full)
summary(nm_envr_recreation)
#                        Estimate     SE     z  P(>|z|)
#(Intercept)                 0.536 0.1751  3.06 2.19e-03
#habitat_simpleRocks         2.397 0.2105 11.39 4.96e-30
#habitat_simpleDwarf_pine   -1.902 0.6081 -3.13 1.77e-03
#habitat_simpleForest       -1.137 0.2186 -5.20 1.98e-07
#elevation_sc               -0.782 0.0823 -9.50 2.09e-21 #0.0000000000000209 #QUESTION!!! WHY WOULD IT BE NEGATIVE SINCE ALPINE SPECIES? PERHAPSH WITHIN HABITAT THEY PREFER LOWER ELEVATIONS? 
#trail_500_sc               -0.420 0.0979 -4.29 1.81e-05 #0.0000181
# Call:
#   pcount(formula = ~log_effort ~ habitat_simple + elevation_sc + 
#            trail_500_sc, data = umf_chamois_count, K = 100)
# Detection (logit-scale):
#   Estimate    SE     z  P(>|z|)
# (Intercept)    -5.07 0.609 -8.33 8.22e-17
# log_effort      1.79 0.313  5.72 1.06e-08
# 
# AIC: 1922.96 
# Number of sites: 104



##CHAMOIS TEMPORAL
chm_data_chamois <- make_chm_data(
  deployments  = deployment_covariates_chm,
  observations = observations_chm %>%
    filter(scientificName == "Rupicapra rupicapra"),
  nBins = 24,
  covs = c("sheep_rate_100tn", "total_trail_m", "habitat_simple"),
  collapse = TRUE) %>%
  mutate(sheep_rate_sc = as.numeric(scale(sheep_rate_100tn)),
         trail_500_sc  = as.numeric(scale(total_trail_m)),
         habitat_simple = case_when(
           locationName == "1RET022" ~ "Open_pasture",
           locationName == "1RET035" ~ "Dwarf_pine",
           TRUE ~ as.character(habitat_simple)) %>% 
           factor(levels = c("Open_pasture", "Rocks",
                             "Dwarf_pine", "Forest")))

length(unique(chm_data_chamois$locationName))   #104

#base models
chm_null_chamois <- mixed_model(
  fixed  = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_chamois)

chm_unimodal_chamois <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_chamois)

chm_bimodal_chamois <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_chamois)

#habitat control
chm_bimodal_hab_chamois <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian) +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_chamois)

#disturbance, habitat controlled
chm_sheep_hab_chamois <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_chamois)

chm_trail_hab_chamois <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_chamois)

chm_trail_sheep_hab_chamois <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    cos(timeRadian) * sheep_rate_sc +
    sin(timeRadian) * sheep_rate_sc +
    cos(2*timeRadian) * sheep_rate_sc +
    sin(2*timeRadian) * sheep_rate_sc +
    habitat_simple,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_chamois)

AIC(chm_null_chamois, chm_unimodal_chamois, chm_bimodal_chamois,
    chm_bimodal_hab_chamois, chm_sheep_hab_chamois, chm_trail_hab_chamois, chm_trail_sheep_hab_chamois)
#                        df      AIC
#chm_null_chamois         2 1520.599
#chm_unimodal_chamois     4 1494.575
#chm_bimodal_chamois      6 1413.104
#chm_bimodal_hab_chamois  9 1394.900 #Δ3.09 worse than best
#chm_sheep_hab_chamois   14 1397.974
#chm_trail_hab_chamois   14 1391.617 #best
#chm_trail_sheep_hab_chamois 19 1394.702 
summary(chm_trail_hab_chamois)
# Fixed effects:
#                                 Estimate Std.Err  z-value    p-value
#(Intercept)                       -8.3084  0.5108 -16.2643    < 1e-04
#cos(timeRadian)                   -0.5155  0.1045  -4.9349    < 1e-04
#trail_500_sc                      -0.3922  0.2760  -1.4212 0.15525433
#sin(timeRadian)                    0.2528  0.0739   3.4231 0.00061903
#cos(2 * timeRadian)               -0.6787  0.0879  -7.7206    < 1e-04
#sin(2 * timeRadian)                0.1186  0.0839   1.4128 0.15770326
#habitat_simpleRocks                2.0428  0.7377   2.7692 0.00561993
#habitat_simpleDwarf_pine          -2.1966  0.9455  -2.3231 0.02017271
#habitat_simpleForest              -0.6774  0.6187  -1.0948 0.27358846
#cos(timeRadian):trail_500_sc      -0.0775  0.1367  -0.5669 0.57078274
#trail_500_sc:sin(timeRadian)       0.2920  0.0973   3.0017 0.00268492
#trail_500_sc:cos(2 * timeRadian)  -0.0907  0.1153  -0.7867 0.43144459
#trail_500_sc:sin(2 * timeRadian)   0.2056  0.1099   1.8702 0.06146244 
summary(chm_sheep_hab_chamois)

#fit_chm comparason
fitchm_trail_habint_chamois <- fit_chm(
  cbind(success, failure) ~ trail_500_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_chamois)

fitchm_sheep_habint_chamois <- fit_chm(
  cbind(success, failure) ~ sheep_rate_sc + habitat_simple,
  type = "bimodal",
  data = chm_data_chamois)

AIC(chm_bimodal_hab_chamois, chm_trail_hab_chamois, chm_sheep_hab_chamois,
    fitchm_trail_habint_chamois, fitchm_sheep_habint_chamois)
# df      AIC
# chm_bimodal_hab_chamois      9 1394.900
# chm_trail_hab_chamois       14 1391.617
# chm_sheep_hab_chamois       14 1397.974
# fitchm_trail_habint_chamois 26 1367.917
# fitchm_sheep_habint_chamois 26 1384.829
summary(fitchm_trail_habint_chamois) 
#fit_chm is reports stronger, chamois do use habitats differently across the day.

##ROE DEER SPATIAL
#rocks sites excluded (0 detections across 14 sites), habitat as Forest vs Non-forest
# site_covs_roe <- site_covs %>%
#   filter(!site_id %in% rocks_sites) %>%
#   mutate(habitat_roe = case_when(
#     habitat_simple == "Forest" ~ "Forest",
#     TRUE ~ "Non-forest") %>% factor(levels = c("Non-forest", "Forest")))
# 
# #trail_500_sc comes through from site_covs, confirm existeance
# summary(site_covs_roe$trail_500_sc)
# table(site_covs_roe$habitat_roe)   # 41 Non-forest, 49 Forest
# 
# site_covs_umf_roe <- site_covs_roe %>%
#   select(habitat_roe, sheep_rate_sc, trail_500_sc, elevation_sc)
# 
# umf_roe_model <- unmarkedFramePCount(
#   y        = as.matrix(roe_deer_matrix_model),
#   siteCovs = as.data.frame(site_covs_umf_roe),
#   obsCovs  = list(log_effort = log_effort_roe_model))

# nm0_roe <- pcount(~ log_effort ~ 1, data = umf_roe_model, K = 100)
# nm_envr_roe <- pcount(~ log_effort ~ habitat_roe + elevation_sc,
#                       data = umf_roe_model, K = 100)
# nm_pastoral_roe <- pcount(~ log_effort ~ sheep_rate_sc,
#                           data = umf_roe_model, K = 100)
# nm_recreation_roe <- pcount(~ log_effort ~ trail_500_sc,
#                             data = umf_roe_model, K = 100)
# nm_envr_pastoral_roe <- pcount(~ log_effort ~ habitat_roe + elevation_sc +
#                                  sheep_rate_sc, data = umf_roe_model, K = 100)
# nm_envr_recreation_roe <- pcount(~ log_effort ~ habitat_roe + elevation_sc +
#                                    trail_500_sc, data = umf_roe_model, K = 100)
# nm_recreation_pastoral_roe <- pcount(~ log_effort ~ trail_500_sc +
#                                        sheep_rate_sc, data = umf_roe_model, K = 100)
# nm_full_roe <- pcount(~ log_effort ~ habitat_roe + elevation_sc +
#                         trail_500_sc + sheep_rate_sc, data = umf_roe_model, K = 100)
# nm_table_roe <- modSel(nm_model_list_roe)
# print(nm_table_roe)
# #                   nPars     AIC delta   AICwt cumltvWt
# #full                    7 1218.39  0.00 9.1e-01     0.91
# #envr+pastoral           6 1224.10  5.71 5.2e-02     0.96
# #envr+recreation         6 1224.84  6.45 3.6e-02     1.00
# #envr                    5 1230.55 12.16 2.1e-03     1.00
# #recreation+pastoral     5 1245.77 27.38 1.0e-06     1.00
# #pastoral                4 1256.87 38.48 4.0e-09     1.00
# #recreation              4 1271.01 52.62 3.4e-12     1.00
# #null                    3 1282.23 63.84 1.2e-14     1.00
# summary(nm_full_roe)
# #                  Estimate     SE        z  P(>|z|)
# #(Intercept)       -0.00157 0.2565 -0.00611 9.95e-01
# #habitat_roeForest  1.21566 0.2512  4.83850 1.31e-06
# #elevation_sc       0.02688 0.0914  0.29403 7.69e-01
# #trail_500_sc      -0.28171 0.1095 -2.57377 1.01e-02
# #sheep_rate_sc     -1.01169 0.5617 -1.80096 7.17e-02

umf_roe <- unmarkedFramePCount(
  y        = as.matrix(roe_deer_matrix), #roe_deer_matrix (104 rows) vs roe_deer_matrix_model (90 rows)
  siteCovs = as.data.frame(site_covs_umf),
  obsCovs  = list(log_effort = log_effort_roe)) #log_Effort_roe (104 sites) vs log_Effort_roe_model (90 sites), will have to sort scripts.

nm0_roe <- pcount(~ log_effort ~ 1, data = umf_roe, K = 100)
nm_envr_roe <- pcount(~ log_effort ~ habitat_simple + elevation_sc,
                      data = umf_roe, K = 100)
nm_pastoral_roe <- pcount(~ log_effort ~ sheep_rate_sc,
                          data = umf_roe, K = 100)
nm_recreation_roe <- pcount(~ log_effort ~ trail_500_sc,
                            data = umf_roe, K = 100)
nm_envr_pastoral_roe <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                 sheep_rate_sc, data = umf_roe, K = 100)
nm_envr_recreation_roe <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                                   trail_500_sc, data = umf_roe, K = 100)
nm_recreation_pastoral_roe <- pcount(~ log_effort ~ trail_500_sc +
                                       sheep_rate_sc, data = umf_roe, K = 100)
nm_full_roe <- pcount(~ log_effort ~ habitat_simple + elevation_sc +
                        trail_500_sc + sheep_rate_sc, data = umf_roe, K = 100)

nm_model_list_roe <- fitList(
  "null"                = nm0_roe,
  "envr"                = nm_envr_roe,
  "pastoral"            = nm_pastoral_roe,
  "recreation"          = nm_recreation_roe,
  "envr+pastoral"       = nm_envr_pastoral_roe,
  "envr+recreation"     = nm_envr_recreation_roe,
  "recreation+pastoral" = nm_recreation_pastoral_roe,
  "full"                = nm_full_roe)

nm_table_roe <- modSel(nm_model_list_roe)
print(nm_table_roe)
#                     nPars     AIC  delta   AICwt cumltvWt
# full                    9 1222.36   0.00 8.8e-01     0.88
# envr+pastoral           8 1227.63   5.27 6.3e-02     0.94
# envr+recreation         8 1228.06   5.71 5.1e-02     1.00
# envr                    7 1232.76  10.40 4.9e-03     1.00
# recreation+pastoral     5 1295.50  73.14 1.2e-16     1.00
# pastoral                4 1304.50  82.14 1.3e-18     1.00
# recreation              4 1313.90  91.54 1.2e-20     1.00
# null                    3 1323.38 101.02 1.0e-22     1.00
summary(nm_full_roe)
# #                        Estimate       SE       z  P(>|z|)
# (Intercept)               -0.0396   0.3376 -0.1171 0.906743
# habitat_simpleRocks      -13.2354 218.4951 -0.0606 0.951697
# habitat_simpleDwarf_pine   0.0739   0.4176  0.1771 0.859452
# habitat_simpleForest       1.2557   0.3401  3.6920 0.000223
# elevation_sc               0.0256   0.0919  0.2782 0.780855
# trail_500_sc              -0.2782   0.1111 -2.5052 0.012237 #sig
# sheep_rate_sc             -1.0039   0.5634 -1.7819 0.074760 #not sig
# Call:
#   pcount(formula = ~log_effort ~ habitat_simple + elevation_sc + 
#            trail_500_sc + sheep_rate_sc, data = umf_roe, K = 100)
# Detection (logit-scale):
#   Estimate    SE     z  P(>|z|)
# (Intercept)   -4.200 0.507 -8.28 1.21e-16
# log_effort     0.871 0.250  3.49 4.89e-04
# 
# AIC: 1222.359 
# Number of sites: 104

###CONCLUSION FOR SWITCH: Roe deer abundance was modelled across all 104 sites.
#No roe deer were detected at any of the 14 rocky sites, so the Rocks
#coefficient is not identifiable (β = −13.24, SE = 218.5); a model excluding
#these sites and collapsing habitat to Forest/Non-forest produced near-identical
#estimates for all other parameters (differences <0.01), confirming that the
#remaining coefficients are unaffected. The full model was clearly best
#supported (ΔAIC = 5.27, weight = 0.88). Roe deer were substantially more
#abundant in forest (β = 1.26, p < 0.001), with no effect of elevation (p =
#0.781) or dwarf pine (p = 0.859). Abundance declined significantly with trail
#density within 500m (β = −0.28, p = 0.012), indicating spatial avoidance of
#recreational infrastructure. Sheep detection rate showed a large negative
#estimate that fell short of significance (β = −1.00, p = 0.075), with a wide
#standard error reflecting co-occurrence at only three camera sites; this should
#be interpreted with caution.

##ROE DEER TEMPORAL

#sheep excluded, roe deer and sheep co-occur at only 3 sites
#rocks sites excluded and habitat as Forest/Non-forest, matching spatial model
chm_data_roe <- make_chm_data(
  deployments  = deployment_covariates_chm %>%
    filter(!locationName %in% rocks_sites),
  observations = observations_chm %>%
    filter(scientificName == "Capreolus capreolus"),
  nBins  = 24,
  covs = c("total_trail_m", "habitat_simple"),
  collapse = TRUE) %>%
  mutate(trail_500_sc = as.numeric(scale(total_trail_m)),
         habitat_simple = case_when(
           locationName == "1RET022" ~ "Open_pasture",
           locationName == "1RET035" ~ "Dwarf_pine",
           TRUE ~ as.character(habitat_simple)),
         habitat_roe = case_when(
           habitat_simple == "Forest" ~ "Forest",
           TRUE ~ "Non-forest") %>% factor(levels = c("Non-forest", "Forest")))

length(unique(chm_data_roe$locationName))   #90
table(chm_data_roe$habitat_roe)

chm_null_roe <- mixed_model(
  fixed  = cbind(success, failure) ~ 1,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_roe)

chm_unimodal_roe <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_roe)

chm_bimodal_roe <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian),
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_roe)

chm_bimodal_hab_roe <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) + sin(timeRadian) +
    cos(2*timeRadian) + sin(2*timeRadian) +
    habitat_roe,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_roe)

chm_trail_hab_roe <- mixed_model(
  fixed  = cbind(success, failure) ~
    cos(timeRadian) * trail_500_sc +
    sin(timeRadian) * trail_500_sc +
    cos(2*timeRadian) * trail_500_sc +
    sin(2*timeRadian) * trail_500_sc +
    habitat_roe,
  random = ~ 1 | locationName,
  family = binomial(),
  data   = chm_data_roe)

AIC(chm_null_roe, chm_unimodal_roe, chm_bimodal_roe,
    chm_bimodal_hab_roe, chm_trail_hab_roe)
#                    df      AIC
#chm_null_roe         2 1314.791
#chm_unimodal_roe     4 1288.863
#chm_bimodal_roe      6 1262.732
#chm_bimodal_hab_roe  7 1240.269
#chm_trail_hab_roe   12 1243.501 
summary(chm_trail_hab_roe)

#fit_chm comparason
fitchm_trail_habint_roe <- fit_chm(
  cbind(success, failure) ~ trail_500_sc + habitat_roe,
  type = "bimodal",
  data = chm_data_roe)

AIC(chm_bimodal_hab_roe,
    chm_trail_hab_roe,
    fitchm_trail_habint_roe)
# df      AIC
# chm_bimodal_hab_roe      7 1240.269
# chm_trail_hab_roe       12 1243.501
# fitchm_trail_habint_roe 16 1247.309
summary(chm_bimodal_hab_roe)
# Call:
#   mixed_model(fixed = cbind(success, failure) ~ cos(timeRadian) + 
#                 sin(timeRadian) + cos(2 * timeRadian) + sin(2 * timeRadian) + 
#                 habitat_roe, random = ~1 | locationName, data = chm_data_roe, 
#               family = binomial())
# Data Descriptives:
#   Number of Observations: 2160
# Number of Groups: 90 
# Model:
#   family: binomial
# link: logit 
# Fit statistics:
#   log.Lik      AIC      BIC
# -613.1345 1240.269 1257.768
# Random effects covariance matrix:
#   StdDev
# (Intercept) 1.241925
# Fixed effects:
#   Estimate Std.Err  z-value   p-value
# (Intercept)          -8.8756  0.3487 -25.4525   < 1e-04
# cos(timeRadian)      -0.6721  0.1218  -5.5186   < 1e-04
# sin(timeRadian)       0.2583  0.0939   2.7503 0.0059545
# cos(2 * timeRadian)  -0.4524  0.1032  -4.3820   < 1e-04
# sin(2 * timeRadian)   0.3246  0.1018   3.1884 0.0014307
# habitat_roeForest     1.8604  0.3787   4.9132   < 1e-04
# Integration:
#   method: adaptive Gauss-Hermite quadrature rule
# quadrature points: 11
# Optimization:
#   method: hybrid EM and quasi-Newton
# converged: TRUE 

#-spatial: forest specialists, trail avoidant. no sheep effect.
#-temporal: strong habitat effect, no trail effect. no sheep inclucded. 

#the full model was best supported, but the sheep coefficient fell short of
#significance with a wide standard error reflecting minimal co-occurrence.roe
#deer and sheep barely co-occur — 3 sites, 5 detections — so this study cannot
#assess whether sheep affect roe deer.

############################### PLOTS TEMPORAL #########################################

library(ggplot2)

##CHAMOIS TEMPORAL - ##chm_trail_hab_chamois AIC 1391.617 #best
newdat_chamois <- expand.grid(
  timeRadian = seq(0, 24, len = 100) *pi/12,
  trail_500_sc = quantile(chm_data_chamois$trail_500_sc, c(0.1, 0.9)),
  habitat_simple = factor("Rocks", levels = levels(chm_data_chamois$habitat_simple)))
#FYI:chamois overlap with sheep far more than roe deer do. 10 of 17 sheep
#sites versus 3, and 82 detections versus 5.

pred_chamois <- effectPlotData(chm_trail_hab_chamois, newdat_chamois,
                               marginal = TRUE) %>% 
  mutate(Trail = ifelse(trail_500_sc == min(trail_500_sc), "Low", "High") %>%
           factor(levels = c("Low", "High")))

#meter values for the legend labels
quantile(chm_data_chamois$total_trail_m, c(0.1, 0.9))

p_chamois <-ggplot(pred_chamois, aes(timeRadian, plogis(pred))) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp),
                  fill = Trail), alpha = 0.2) +
  geom_line(aes(colour = Trail), linewidth = 1) +
  labs(x = "Time of day (hour)",
       y = "Predicted activity (probability)",
       title = "Chamois") +
  scale_x_continuous(breaks = seq(0, 2*pi, len = 5),
                     labels = seq(0, 24, len = 5)) + theme_minimal() +
  theme(legend.position = "top")
p_chamois
#Why chamois held at Rocks? Rocks because chamois are rock specialists (84% of
#detections), so it's their realistic context. The levels =  part matters: the
#factor must have all four levels defined even though you're using one, or the
#model won't recognise it. Since habitat is additive in your model, this choice
#only moves the curves vertically. The gap between the low and high trail
#curves, the thing the figure shows, is identical whichever level I pick. 

#How many sites have any trail within 500m, and how many chamois detections at
#those sites? #26 out of 58 trail sites, with 208 detections.
site_covs %>% mutate(chamois_det = rowSums(chamois_count_matrix, na.rm = TRUE),
                     has_trail   = total_trail_m > 0) %>% group_by(has_trail) %>%
  summarise(n_sites = n(), sites_with_chamois = sum(chamois_det > 0),
            total_detections = sum(chamois_det))


##BEAR TEMPORAL - ##chm_sheep_hab_bear AIC 1221.387 #best
newdat_bear <- expand.grid(
  timeRadian = seq(0, 24, len = 100) *pi/12,
  sheep_rate_sc = quantile(chm_data_bear$sheep_rate_sc, c(0.1, 0.9)),
  habitat_simple = factor("Open_pasture",
                          levels = levels(chm_data_bear$habitat_simple)))

pred_bear <- effectPlotData(chm_sheep_hab_bear,
                            newdat_bear, marginal = TRUE) %>%
  mutate(Sheep = ifelse(sheep_rate_sc == min(sheep_rate_sc), "Low", "High") %>%
           factor(levels = c("Low", "High")))

quantile(chm_data_bear$sheep_rate_100tn, c(0.1, 0.9))

p_bear <- ggplot(pred_bear, aes(timeRadian, plogis(pred))) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp),
                  fill = Sheep), alpha = 0.2) +
  geom_line(aes(colour = Sheep), linewidth = 1) +
  labs(x = "Time of day (hour)",
       y = "Predicted activity (probability)",
       title = "Brown bear") + scale_x_continuous(breaks = seq(0, 2*pi, len = 5),
                                                  labels = seq(0, 24, len = 5)) + theme_minimal() +
  theme(legend.position = "top")
p_bear

#Why bears held at open pasture? 
#How many sites have sheep detecs and how many bear detections at those sites? 
site_covs %>% mutate(bear_det  = rowSums(bear_matrix, na.rm = TRUE),
                     has_sheep = sheep_rate_100tn > 0) %>%group_by(has_sheep) %>%
  summarise(n_sites = n(), sites_with_bear = sum(bear_det > 0),
            total_detections = sum(bear_det))
#Bears are most abundant there, in your bear spatial model, Open_pasture is the
#reference level (intercept 1.217), and Rocks (−1.12, p = 0.0017) and Forest
#(−0.86, p = 0.0012) are both significantly lower. Dwarf pine is
#indistinguishable. And sheep mostly in open pasture.


##REDDEER TEMPORAL - ##chm_sheep_hab_rd AIC 1517.720 #best
newdat_rd <- expand.grid(
  timeRadian = seq(0, 24, len = 100) *pi/12,
  sheep_rate_sc = quantile(chm_data_red_deer$sheep_rate_sc, c(0.1, 0.9)),
  habitat_simple = factor("Open_pasture", levels = levels(chm_data_red_deer$habitat_simple)))

pred_rd <- effectPlotData(chm_sheep_hab_rd, newdat_rd, marginal = TRUE) %>%
  mutate(Sheep = ifelse(sheep_rate_sc == min(sheep_rate_sc), "Low", "High") %>%
           factor(levels = c("Low", "High")))

quantile(chm_data_red_deer$sheep_rate_100tn, c(0.1, 0.9))

p_rd <- ggplot(pred_rd, aes(timeRadian, plogis(pred))) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp),
                  fill = Sheep), alpha = 0.2) +
  geom_line(aes(colour = Sheep), linewidth = 1) +
  labs(x = "Time of day (hour)",
       y = "Predicted activity (probability)",
       title = "Red deer") +
  scale_x_continuous(breaks = seq(0, 2*pi, len = 5), labels = seq(0, 24, len = 5)) +
  theme_minimal() + theme(legend.position = "top")
p_rd
#Red deer delayed their evening activity peak at high-sheep sites.

#Why reddeer held at open pasture? 
#How many sites have sheep detecs and how many rd detections at those sites? 
site_covs %>% mutate(rd_det = rowSums(red_deer_matrix, na.rm = TRUE),
                     has_sheep = sheep_rate_100tn > 0) %>% group_by(has_sheep) %>%
  summarise(n_sites = n(), sites_with_rd = sum(rd_det > 0),
            total_detections = sum(rd_det))
#Red deer × sheep: 8 of 17 sheep sites have red deer, with 55 detections.
#Thinner than bear (13 sites, 36 detections) in site count but more detections.
#Enough to fit on, and nothing like roe deer's 3 sites. 

##FOX TEMPORAL - ##chm_sheep_hab_fox AIC 1178.928 #best
newdat_fox <- expand.grid(
  timeRadian = seq(0, 24, len = 100) *pi/12,
  sheep_rate_sc = quantile(chm_data_fox$sheep_rate_sc, c(0.1, 0.9)),
  habitat_simple = factor("Open_pasture", levels = levels(chm_data_fox$habitat_simple)))

pred_fox <- effectPlotData(chm_sheep_hab_fox, newdat_fox,
                           marginal = TRUE) %>% 
  mutate(Sheep = ifelse(sheep_rate_sc == min(sheep_rate_sc), "Low", "High") %>%
           factor(levels = c("Low", "High")))

quantile(chm_data_fox$sheep_rate_100tn, c(0.1, 0.9))

p_fox <- ggplot(pred_fox, aes(timeRadian, plogis(pred))) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp), fill = Sheep), alpha = 0.2) +
  geom_line(aes(colour = Sheep), linewidth = 1) +
  labs(x = "Time of day (hour)",
       y = "Predicted activity (probability)",
       title = "Red fox") + scale_x_continuous(breaks = seq(0, 2*pi, len = 5),
                                               labels = seq(0, 24, len = 5)) + theme_minimal() + theme(legend.position = "top")
p_fox

##Why fox held at open pasture? both abudant at open pasture.
#How many sites have sheep detecs and how many fox detections at those sites?
#Very abudant at sheep sites.
site_covs %>% mutate(fox_det = rowSums(red_fox_matrix, na.rm = TRUE),
                     has_sheep = sheep_rate_100tn > 0) %>% group_by(has_sheep) %>%
  summarise(n_sites = n(), sites_with_fox = sum(fox_det > 0),
            total_detections = sum(fox_det))

##ROE DEER TEMPORAL - ##
#fixed at which habitat?
#How many sites have sheep detecs and how many roe deer detections at those sites?
site_covs %>% mutate(roe_det = rowSums(roe_deer_matrix, na.rm = TRUE),
                     has_sheep = sheep_rate_100tn > 0) %>%
  group_by(has_sheep) %>% summarise(n_sites = n(), sites_with_roe = sum(roe_det > 0),
                                    total_detections = sum(roe_det)) #co-occured at 3 sites, 5 roe deer detections.
#roe deer detections mostly in 

newdat_roe <- data.frame(
  timeRadian  = seq(0, 24, len = 100) *pi/12,
  habitat_roe = factor("Forest", levels = levels(chm_data_roe$habitat_roe)))

pred_roe <- effectPlotData(chm_bimodal_hab_roe, newdat_roe, marginal = TRUE)

p_roe <- ggplot(pred_roe, aes(timeRadian, plogis(pred))) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp)),fill = "grey60", alpha = 0.3) +
  geom_line(linewidth = 1) + labs(x = "Time of day (hour)",
                                  y = "Predicted activity (probability)",
                                  title = "Roe deer") +
  scale_x_continuous(breaks = seq(0, 2*pi, len = 5),labels = seq(0, 24, len = 5)) +
  theme_minimal()
p_roe


#SHEEP
newdat_sheep <- data.frame(timeRadian = seq(0, 24, len = 100) *pi/12)
pred_sheep <- effectPlotData(chm_bimodal_sheep, newdat_sheep, marginal = TRUE)
max(plogis(pred_sheep$pred))   #compare to ~0.003 for bear

summary(chm_bimodal_sheep)
range(plogis(pred_sheep$pred))

############################### PLOTS SPATIAL #########################################

##CHAMOIS SPATIAL - ###envr+recreation AIC 1922.96 #best
#Fixed habitat: Rocks (same as temporal)
trail_seq <- seq(min(site_covs$trail_500_sc), max(site_covs$trail_500_sc),
                 length.out = 100)

newdat_chamois_sp <- data.frame(trail_500_sc = trail_seq,elevation_sc = 0,
                                habitat_simple = factor("Rocks", levels = levels(site_covs$habitat_simple)))

pred_chamois_sp <- unmarked::predict(nm_envr_recreation,
                                     type = "state", newdata = newdat_chamois_sp,
                                     appendData = TRUE) %>%
  mutate(trail_m = trail_500_sc * sd(site_covs$total_trail_m) + mean(site_covs$total_trail_m))

p_chamois_sp <- ggplot(pred_chamois_sp, aes(trail_m, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60", alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(x = "Trail density within 500 m (m)",
       y = "Predicted abundance",
       title = "Chamois") + theme_minimal()
p_chamois_sp


##BEAR SPATIAL - #envr+recreation AIC 1187.57 #best
#Fixed habitat: Open pasture (same as temporal)
newdat_bear_sp <- data.frame(trail_500_sc = trail_seq, elevation_sc = 0,
                             habitat_simple = factor("Open_pasture", levels = levels(site_covs$habitat_simple)))

pred_bear_sp <- unmarked::predict(nm_envr_recreation_bear,
                                  type = "state", newdata = newdat_bear_sp, appendData = TRUE) %>%
  mutate(trail_m = trail_500_sc * sd(site_covs$total_trail_m) +
           mean(site_covs$total_trail_m))

p_bear_sp <- ggplot(pred_bear_sp, aes(trail_m, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60", alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(x = "Trail density within 500 m (m)",
       y = "Predicted abundance",
       title = "Brown bear") + theme_minimal()
p_bear_sp

##REDDEER SPATIAL - ###envr+recreation AIC 1649.58 #best
#Fixed habitat: Forest (not same as temporal)
newdat_rd_sp <- data.frame(trail_500_sc = trail_seq,elevation_sc = 0,
                           habitat_simple = factor("Forest", levels = levels(site_covs$habitat_simple)))

pred_rd_sp <- unmarked::predict(nm_envr_recreation_rd, type = "state",
                                newdata = newdat_rd_sp, appendData = TRUE) %>%
  mutate(trail_m = trail_500_sc * sd(site_covs$total_trail_m) +
           mean(site_covs$total_trail_m))

p_rd_sp <- ggplot(pred_rd_sp, aes(trail_m, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60", alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(x = "Trail density within 500 m (m)",
       y = "Predicted abundance",
       title = "Red deer") + theme_minimal()
p_rd_sp

##ROE SPATIAL - ###full AIC 1222.36 #best #but sheep not sig, so plotting trail.
#Fixed habitat: Forest (same as temporal)
newdat_roe_sp <- data.frame(trail_500_sc   = trail_seq, elevation_sc   = 0,
                            sheep_rate_sc  = 0,
                            habitat_simple = factor("Forest", levels = levels(site_covs$habitat_simple)))

pred_roe_sp <- unmarked::predict(nm_full_roe, type = "state", newdata = newdat_roe_sp,
                                 appendData = TRUE) %>% mutate(trail_m = trail_500_sc * sd(site_covs$total_trail_m) +
                                                                 mean(site_covs$total_trail_m))

p_roe_sp <- ggplot(pred_roe_sp, aes(trail_m, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60", alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(x = "Trail density within 500 m (m)",
       y = "Predicted abundance",
       title = "Roe deer") + theme_minimal()
p_roe_sp

##FOX SPATIAL - ###envr+pastoral AIC 1067.98 #best, p-value for sheep not sig
#Fixed habitat: open pasture (not same as temporal)
newdat_fox_sp <- data.frame(trail_500_sc   = trail_seq,elevation_sc   = 0,
                            habitat_simple = factor("Open_pasture", levels = levels(site_covs$habitat_simple)))

pred_fox_sp <- unmarked::predict(nm_envr_recreation_red_fox, type = "state",
                                 newdata = newdat_fox_sp,
                                 appendData = TRUE) %>%
  mutate(trail_m = trail_500_sc * sd(site_covs$total_trail_m) +
           mean(site_covs$total_trail_m))

p_fox_sp <- ggplot(pred_fox_sp, aes(trail_m, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60", alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(x = "Trail density within 500 m (m)",
       y = "Predicted abundance",
       title = "Red fox") + theme_minimal()
p_fox_sp
