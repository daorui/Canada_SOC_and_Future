library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(h2o)

model_dir <- 'models/mineral_bd/'
map_dir <- 'maps/'
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

## step1 training data

mineral_sf <- readRDS('data/processed/mineral_bd_splined.rds')
soc_mineral <- readRDS('data/processed/mineral_soc_splined.rds')
mineral_bd_long <- mineral_sf %>% pivot_longer(cols = ends_with('cm'), names_to = 'depth', values_to = 'bd') %>% mutate(Peat = 0) %>% filter(bd > 0, complete.cases(.)) %>% mutate(ID0 = paste0(id, lon, lat, depth))
mineral_soc_long <- soc_mineral %>% pivot_longer(cols = ends_with('cm'), names_to = 'depth', values_to = 'SOC') %>% mutate(Peat = 0) %>% filter(SOC > 0, complete.cases(.)) %>% mutate(ID0 = paste0(id, lon, lat, depth))
mineral_bd_dataset_long <- left_join(mineral_bd_long, mineral_soc_long, by = 'ID0')

t11 <- mineral_bd_dataset_long %>% select(c(3:6, 14:15))
names(t11) <- c('lon', 'lat', 'depth', 'bd', 'SOC', 'Peat')
t111 <- t11 %>% mutate(site = paste0(lon, '_', lat)) %>% group_by(site) %>% filter(row_number() != n()) %>% filter(n() > 3) %>% ungroup()

## step2 covariate extraction and preprocessing

bd_ca_sf0 <- st_as_sf(t111, coords = c('lon', 'lat'), crs = st_crs(4326))
ca_covas <- rast('data/covariates/covariates_bd.tif')
reclass_table <- matrix(c(
  0, 5, 5,
  5, 30, 30,
  30, 45, 45,
  45, 100, 100
), ncol = 3, byrow = TRUE)
peat_reclassified22 <- classify(ca_covas$Peat, reclass_table, include.lowest = T, right = T)
ca_covas[[62]] <- peat_reclassified22
ttt <- terra::extract(ca_covas, vect(bd_ca_sf0))
bd_ca_sf1 <- cbind(bd_ca_sf0[, 1:3], ttt)
bd_ca_sf1$Pop_den28 <- as.integer(bd_ca_sf1$Pop_den28)
bd_ca_sf1$Peat <- factor(bd_ca_sf1$Peat, levels = c(5, 30, 60, 100))

names(bd_ca_sf1)[1] <- 'Depth'
d10 <-
  bd_ca_sf1 %>% select(c(1:3, 5:66)) %>% mutate(
    Depth = case_when(
      Depth == '0-5 cm' ~ 2.5,
      Depth == '5-15 cm' ~ 10,
      Depth == '5-10 cm' ~ 7.5,
      Depth == '10-15 cm' ~ 12.5,
      Depth == '15-20 cm' ~ 17.5,
      Depth == '15-30 cm' ~ 22.5,
      Depth == '20-30 cm' ~ 25,
      Depth == '30-40 cm' ~ 35,
      Depth == '30-45 cm' ~ 37.5,
      Depth == '40-50 cm' ~ 45,
      Depth == '50-60 cm' ~ 55,
      Depth == '45-60 cm' ~ 52.5,
      Depth == '60-70 cm' ~ 65,
      Depth == '60-80 cm' ~ 70,
      Depth == '70-80 cm' ~ 75,
      Depth == '80-90 cm' ~ 85,
      Depth == '90-100 cm' ~ 95,
      Depth == '80-100 cm' ~ 90,
      Depth == '100-120 cm' ~ 110,
      Depth == '120-150 cm' ~ 135,
      Depth == '150-180 cm' ~ 165,
      Depth == '180-210 cm' ~ 195,
      Depth == '210-240 cm' ~ 225,
      Depth == '240-270 cm' ~ 255,
      Depth == '270-300 cm' ~ 285,
    )
  ) %>%
  mutate(
    NDVI1 = NDVI1 / 10000,
    EVI2 = EVI2 / 10000,
    NPP3 = NPP3 / 10000,
    GPP4 = GPP4 / 1000,
    PET5 = log(PET5),
    AI6 = AI6 / 10000,
    BIO4 = BIO4 / 100,
    Solar26 = Solar26 / 10000,
    So_seasonal27 = So_seasonal27 / 10000,
    Vrm30 = Vrm30 * 10^4,
    Tcurv31 = Tcurv31 * 10^4,
    Pcurv32 = Pcurv32 * 10^4,
    Ref_band1 = Ref_band1 / 10000,
    ref_band2 = ref_band2 / 10000,
    ref_band3 = ref_band3 / 10000,
    ref_band4 = ref_band4 / 10000,
    ref_band5 = ref_band5 / 10000,
    ref_band6 = ref_band6 / 10000,
    ref_band7 = ref_band7 / 10000,
    Maj48 = factor(Maj48, levels = 1:10),
    lith.8 = factor(lith.8, levels = 1:7),
    lith.24 = factor(lith.24, levels = 1:15),
    Peat = factor(Peat, levels = c(5, 30, 60, 100)),
    lon = lon / 10^6,
    lat = lat / 10^6
  ) %>% na.omit()
d10[, paste0('BIO', 12:19)] <- log(st_drop_geometry(d10[, paste0('BIO', 12:19)]))
d11 <- d10
d11$Pop_den28 <- as.integer(d11$Pop_den28)
d11 <- d11 %>% filter(Depth < 180)

## step3 model building

h2o.init(nthreads = 32, max_mem_size = '100G')
d1_hdf <- as.h2o(x = sf::st_drop_geometry(d11))
set.seed(23)
h_auto <- h2o.automl(x = setdiff(names(d1_hdf), c('NDVI1', 'EVI2', 'NPP3', 'bd', 'lon', 'lat', 'Peat')), y = 'bd',
                     training_frame = d1_hdf,
                     exclude_algos = c('GLM', 'DeepLearning'),
                     max_models = NULL,
                     nfolds = 10,
                     stopping_metric = 'RMSE',
                     keep_cross_validation_predictions = TRUE,
                     seed = 927)
lb <- as.data.frame(h2o.get_leaderboard(object = h_auto, extra_columns = 'ALL'))
saveRDS(lb, file.path(model_dir, 'leaderboard_mineral_bd.rds'))
stack_m <- h2o.get_best_model(h_auto, algorithm = 'StackedEnsemble', criterion = 'RMSE')
h2o.saveModel(stack_m, model_dir, force = TRUE)

## step4 map prediction

ca_covas <- rast('data/covariates/covariates_bd.tif')
peat_ca <- rast('data/covariates/peat_fraction.tif')
peat_ca <- ifel(is.na(peat_ca), 0, peat_ca)
reclass_table <- matrix(c(
  0, 5, 5,
  5, 30, 30,
  30, 60, 60,
  60, 100, 100
), ncol = 3, byrow = TRUE)
peat_reclassified <- classify(peat_ca, reclass_table, include.lowest = T, right = T)
ca_covas[['Peat']] <- peat_reclassified

dep00 <- c(0, 5, 15, 30, 60, 100)
t2 <- rast(file.path(map_dir, paste0('mineral_soc_', dep00, 'cm.tif')))

ca_covas$Depth <- 0
data0 <- as.data.frame(as.matrix(ca_covas))
ca_base_new <- data0 %>%
  mutate_all(list(~ replace(., . == c(-9999), NA))) %>%
  mutate(
    NDVI1 = NDVI1 / 10000,
    EVI2 = EVI2 / 10000,
    NPP3 = NPP3 / 10000,
    GPP4 = GPP4 / 1000,
    PET5 = log(PET5),
    AI6 = AI6 / 10000,
    BIO4 = BIO4 / 100,
    BIO12 = ifelse(BIO12 == 0, -1, log(BIO12)),
    BIO13 = ifelse(BIO13 == 0, -1, log(BIO13)),
    BIO14 = ifelse(BIO14 == 0, -1, log(BIO14)),
    BIO15 = ifelse(BIO15 == 0, -1, log(BIO15)),
    BIO16 = ifelse(BIO16 == 0, -1, log(BIO16)),
    BIO17 = ifelse(BIO17 == 0, -1, log(BIO17)),
    BIO18 = ifelse(BIO18 == 0, -1, log(BIO18)),
    BIO19 = ifelse(BIO19 == 0, -1, log(BIO19)),
    Solar26 = Solar26 / 10000,
    So_seasonal27 = So_seasonal27 / 10000,
    Vrm30 = Vrm30 * 10^4,
    Tcurv31 = Tcurv31 * 10^4,
    Pcurv32 = Pcurv32 * 10^4,
    Ref_band1 = Ref_band1 / 10000,
    ref_band2 = ref_band2 / 10000,
    ref_band3 = ref_band3 / 10000,
    ref_band4 = ref_band4 / 10000,
    ref_band5 = ref_band5 / 10000,
    ref_band6 = ref_band6 / 10000,
    ref_band7 = ref_band7 / 10000,
    Maj48 = factor(Maj48, levels = 1:10),
    lith.8 = factor(lith.8, levels = 1:7),
    lith.24 = factor(lith.24, levels = 1:15),
    Peat = factor(Peat, levels = c(5, 30, 60, 100)),
    lon = lon / 10^6,
    lat = lat / 10^6
  )
rm(data0); gc()
ind <- rowSums(is.na(ca_base_new)) == 0
pred_df <- ca_base_new[ind, ] %>% mutate(across(where(is.numeric), ~ round(., 4)))
SOC_pred <- as.data.frame(as.matrix(t2))[ind, ]
names(SOC_pred) <- dep00

template <- rast(ca_covas[[1]])
for (i in dep00) {
  pred_df$Depth <- i
  pred_df$SOC <- SOC_pred[, as.character(i)]
  test <- as.h2o(x = pred_df)
  p <- as.data.frame(h2o.predict(stack_m, test))$predict
  p <- pmin(p, 2.5)
  v <- rep(NA_real_, ncell(template)); v[which(ind)] <- p
  r <- setValues(template, v); names(r) <- 'Stack'
  writeRaster(r, file.path(map_dir, paste0('mineral_bd_', i, 'cm.tif')), overwrite = TRUE)
}
h2o.shutdown(prompt = FALSE)
