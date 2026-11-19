library(dplyr)
library(sf)
library(tidyr)
library(aqp)
library(ithir)

raw_dir <- 'data/raw/'
out_dir <- 'data/processed/'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## step1 data loading and merging

use_layer3 <- read.csv(paste0(raw_dir, 'wosis2019_soc_layers.csv'))
use_layer3[90448, 'orgc_value_avg'] <- 188.6
ca_oringal <- use_layer3 %>% filter(country_name == 'Canada') %>% st_as_sf(coords = c('longitude', 'latitude'))

bd <- read.csv(paste0(raw_dir, 'wosis2019_bd_layers_canada.csv'))

soc0 <- ca_oringal[, c(2, 6:8)] %>% mutate(profile_id = as.integer(profile_id))
bd0 <- bd[, c(1, 7, 8, 15)]
soc_bd <- left_join(soc0, bd0, by = c('profile_id', 'upper_depth', 'lower_depth'), multiple = 'first')

peat_use <- readRDS(paste0(raw_dir, 'canadian_peat_profile_database.rds'))
peat_use1 <- peat_use[, c(1:5, 17, 34, 35)] %>% mutate(C_ORG_PCT = C_ORG_PCT * 10)
peat_use1$SAMP_THICK[which(is.na(peat_use1$SAMP_THICK))] <- 3
peat_use2 <- peat_use1 %>% mutate(lower_depth = UPPER_SAMP_DEPTH + SAMP_THICK) %>% filter(complete.cases(.))

wosis <- soc_bd %>% mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2], source = 'wosis') %>% st_drop_geometry()
sasha_peat <- peat_use2 %>% select(-c(SAMPLE_NO, SAMP_THICK)) %>%
  mutate(profile_id = CORE_ID, upper_depth = UPPER_SAMP_DEPTH, lower_depth = lower_depth,
         orgc_value_avg = C_ORG_PCT, bd_use = BULK_DENSITY, lon = LONGITUDE, lat = LATITUDE, source = 'sasha_peat') %>% select(8, 9, 7, 10:14)

soc_bd_relation_all <- rbind(wosis, sasha_peat)

## step2 peat and mineral profile classification

soc_use_r <- soc_bd_relation_all %>% pivot_longer(cols = c('upper_depth', 'lower_depth'), values_to = 'depth')

wosis_peat_id <- soc_use_r %>% filter(orgc_value_avg > 175, source == 'wosis') %>%
  group_by(profile_id) %>%
  reframe(n = n(), depthmax = max(depth)) %>%
  filter(n >= 3, depthmax >= 40) %>%
  pull(profile_id)
wosis_peat_id <- c(wosis_peat_id, 230535)
wosis_peat <- soc_use_r %>% filter(profile_id %in% wosis_peat_id)
soc_mineral <- soc_use_r %>% filter(!profile_id %in% wosis_peat_id, source == 'wosis')
soc_peat <- rbind(wosis_peat, soc_use_r %>% filter(source == 'sasha_peat'))

## step3 removal of organic surface horizons from mineral profiles

forest_floor <-
  soc_bd_relation_all %>%
  filter(profile_id %in% unique(soc_mineral$profile_id)) %>%
  group_by(profile_id) %>%
  reframe(soc_max = max(head(orgc_value_avg, 3)), soc_min = min(head(orgc_value_avg, 3)),
          ratio = max(head(orgc_value_avg, 3)) / min(head(orgc_value_avg, 3)), n = n()) %>%
  filter(soc_max >= 100, ratio > 4, n >= 3)
forest_floor_id <- unique(forest_floor$profile_id)

forest_floor0 <-
  soc_bd_relation_all %>%
  filter(profile_id %in% setdiff(unique(soc_mineral$profile_id), forest_floor_id)) %>% filter(lower_depth <= 35) %>%
  group_by(profile_id) %>%
  reframe(soc_mean = mean(orgc_value_avg), soc_sd = sd(orgc_value_avg)) %>%
  filter((soc_mean > 100 & soc_sd > 150) | soc_mean > 400)
forest_floor_id0 <- unique(forest_floor0$profile_id)

forest_floor_id2 <-
  soc_bd_relation_all %>%
  filter(profile_id %in% setdiff(unique(soc_mineral$profile_id), c(forest_floor_id0, forest_floor_id))) %>%
  group_by(profile_id) %>%
  reframe(SOC_min = min(head(orgc_value_avg, 3)), SOC_mean = mean(head(orgc_value_avg, 3))) %>% filter(SOC_mean > 200) %>% pull(profile_id)

pure_floor_depth <-
  soc_bd_relation_all %>%
  filter(profile_id %in% c(forest_floor_id0, forest_floor_id, forest_floor_id2)) %>% filter(lower_depth <= 36) %>%
  filter(orgc_value_avg >= 100) %>% group_by(profile_id) %>% reframe(lower_max = max(lower_depth))

mineral_floor_depth <-
  soc_bd_relation_all %>%
  filter(profile_id %in% c(forest_floor_id0, forest_floor_id, forest_floor_id2)) %>% left_join(pure_floor_depth) %>%
  filter(lower_depth > lower_max) %>% mutate(upper_depth1 = upper_depth - lower_max,
                                             lower_depth1 = lower_depth - lower_max) %>%
  filter(upper_depth1 >= 0) %>% mutate(upper_depth = upper_depth1, lower_depth = lower_depth1) %>% select(1:8)

pure_mineral_id <- setdiff(unique(soc_mineral$profile_id), c(forest_floor_id0, forest_floor_id, forest_floor_id2))
pure_mineral <- soc_bd_relation_all %>% filter(profile_id %in% pure_mineral_id)

## step4 depth harmonisation of mineral SOC

mineral_dataset0 <- rbind(pure_mineral, mineral_floor_depth[, 1:8]) %>% filter(lower_depth > 0, orgc_value_avg <= 580)
mineral_dataset1 <- mineral_dataset0
depths(mineral_dataset1) <- profile_id ~ upper_depth + lower_depth
use_bd0 <- ithir::ea_spline(mineral_dataset1, var.name = 'orgc_value_avg', lam = 0.1, d = c(0, 5, 15, 30, 45, 60, 80, 100))
mineral_xy <- mineral_dataset0[, c(1, 6, 7)] %>% distinct()
mineral_sf <- left_join(use_bd0[[1]], mineral_xy, by = c('id' = 'profile_id'))
saveRDS(mineral_sf, paste0(out_dir, 'mineral_soc_splined.rds'))

## step5 screening and depth harmonisation of mineral bulk density

mineral_dataset0 <- rbind(pure_mineral, mineral_floor_depth[, 1:8]) %>% filter(lower_depth > 0, orgc_value_avg <= 580)
mineral_dataset0 <- mineral_dataset0 %>% filter(bd_use > 0, !is.na(bd_use))

d <- mineral_dataset0
d_sus <- subset(d, orgc_value_avg < 50 & bd_use < 1)
m_exp <- nls(
  bd_use ~ a * exp(-b * orgc_value_avg) + c,
  data = d_sus,
  start = list(a = 1.1, b = 0.08, c = 0.05),
  control = nls.control(maxiter = 200)
)
d$bd_fitline <- predict(m_exp, newdata = d)
d$resid_line <- d$bd_use - d$bd_fitline
d$flag_formula_bd <- with(
  d,
  orgc_value_avg < 70 &
    bd_use < 1.2 &
    abs(resid_line) < 0.2
)
d4 <- d %>% filter(!flag_formula_bd, bd_use > 0.6) %>% filter(!bd_use < 0.95 | !orgc_value_avg < 8)

mineral_dataset1 <- d4 %>% select(1:8)
depths(mineral_dataset1) <- profile_id ~ upper_depth + lower_depth
use_bd0 <- ithir::ea_spline(mineral_dataset1, var.name = 'bd_use', lam = 0.1, d = c(0, 5, 15, 30, 45, 60, 80, 100, 120))
mineral_xy <- mineral_dataset0[, c(1, 6, 7)] %>% distinct()
mineral_sf <- left_join(use_bd0[[1]], mineral_xy, by = c('id' = 'profile_id'))
saveRDS(mineral_sf, paste0(out_dir, 'mineral_bd_splined.rds'))

## step6 depth harmonisation of peat SOC and bulk density

soc_peat0 <- soc_bd_relation_all %>% filter(profile_id %in% unique(soc_peat$profile_id))
soc_peat1 <- soc_peat0
depths(soc_peat1) <- profile_id ~ upper_depth + lower_depth
use_bd0 <- ithir::ea_spline(soc_peat1, var.name = 'orgc_value_avg', lam = 0.1, d = c(0, 5, 10, 15, seq(20, 100, 10), seq(120, 300, 30)))
peat_xy <- soc_peat0[, c(1, 6, 7)] %>% distinct()
peat_soc_sf <- left_join(use_bd0[[1]], peat_xy, by = c('id' = 'profile_id'))

soc_peat0 <- soc_bd_relation_all %>% filter(profile_id %in% unique(soc_peat$profile_id))
soc_peat0 <- soc_peat0 %>% filter(bd_use > 0, !is.na(bd_use)) %>%
  mutate(bd_use = ifelse(orgc_value_avg > 450 & bd_use > 0.15, 0.061494 + 1.229292 * exp(c(-0.007661) * orgc_value_avg), bd_use))
soc_peat1 <- soc_peat0
depths(soc_peat1) <- profile_id ~ upper_depth + lower_depth
use_bd0 <- ithir::ea_spline(soc_peat1, var.name = 'bd_use', lam = 0.1, d = c(0, 5, 10, 15, seq(20, 100, 10), seq(120, 300, 30)))
peat_xy <- soc_peat0[, c(1, 6, 7)] %>% distinct()
peat_bd_sf <- left_join(use_bd0[[1]], peat_xy, by = c('id' = 'profile_id'))

## step7 peat carbon density in the upper metre

peat_soc_1m <- peat_soc_sf %>% filter(`soil depth` >= 100) %>% select(1:13)
peat_bd_1m <- peat_bd_sf %>% filter(`soil depth` >= 100) %>% select(1:13)

peat_socd <-
  left_join(peat_soc_1m, peat_bd_1m, by = 'id') %>%
  mutate(`0-5 cm.y` = ifelse(is.na(`0-5 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `0-5 cm.x`), `0-5 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`5-10 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `5-10 cm.x`), `5-10 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`10-15 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `10-15 cm.x`), `10-15 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`15-20 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `15-20 cm.x`), `15-20 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`20-30 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `20-30 cm.x`), `20-30 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`30-40 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `30-40 cm.x`), `30-40 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`40-50 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `40-50 cm.x`), `40-50 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`50-60 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `50-60 cm.x`), `50-60 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`60-70 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `60-70 cm.x`), `60-70 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`70-80 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `70-80 cm.x`), `70-80 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`80-90 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `80-90 cm.x`), `80-90 cm.y`),
         `0-5 cm.y` = ifelse(is.na(`90-100 cm.y`), 0.061494 + 1.229292 * exp(c(-0.007661) * `90-100 cm.x`), `90-100 cm.y`)) %>%
  mutate(`socd1m` = (
    pick(2) * pick(14) * 5 +
      pick(3) * pick(15) * 5 +
      pick(4) * pick(16) * 5 +
      pick(5) * pick(17) * 5 +
      pick(6) * pick(18) * 10 +
      pick(7) * pick(19) * 10 +
      pick(8) * pick(20) * 10 +
      pick(9) * pick(21) * 10 +
      pick(10) * pick(22) * 10 +
      pick(11) * pick(23) * 10 +
      pick(12) * pick(24) * 10 +
      pick(13) * pick(25) * 10
  ) / 100) %>% select(c(1, 26))

peat_socd_use <- data.frame(id = peat_socd$id, socd1m = peat_socd$socd1m$`0-5 cm.x`)
peat_socd_xy <- peat_soc_sf[, c(1, 22, 23)] %>% distinct()
peat_socd_1m <- left_join(peat_socd_use, peat_socd_xy, by = 'id')
saveRDS(peat_socd_1m, paste0(out_dir, 'peat_socd_1m.rds'))
