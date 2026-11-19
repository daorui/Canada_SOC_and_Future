library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(blockCV)
library(h2o)

model_dir <- 'models/peat_socd/'
map_dir <- 'maps/'
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

## step1 training data and covariate extraction

ca_covas <- rast('data/covariates/covariates_peat.tif')
soc_peat <- readRDS('data/processed/peat_socd_1m.rds')
soc_sf <- st_as_sf(soc_peat, coords = c('lon', 'lat'), crs = 4326)
soc_cov <- terra::extract(ca_covas, soc_sf)
soc_u <- cbind(soc_sf, soc_cov)

## step2 spatial cross-validation folds

spatialid <- soc_u %>% group_by(id) %>% reframe(geometry = geometry[1]) %>% st_as_sf()
ca_data_u <- spatialid
set.seed(123)
sb <- cv_spatial(x = ca_data_u, hexagon = TRUE, size = 80000, k = 10, selection = 'random', iteration = 100, biomod2 = FALSE, plot = FALSE, progress = FALSE)
ca_data_u$fold_col <- sb$folds_ids
soc_u2 <- left_join(soc_u, st_drop_geometry(ca_data_u), by = 'id')

## step3 covariate preprocessing

d1 <- soc_u2
covs0 <- c('NDVI1','EVI2','NPP3','GPP4','PET5','AI6','BIO1','BIO10','BIO11','BIO12','BIO13','BIO14','BIO15','BIO16',
           'BIO17','BIO18','BIO19','BIO2','BIO3','BIO4','BIO5','BIO6','BIO7','BIO8','BIO9','Solar26','So_seasonal27',
           'Pop_den28','HFP29','Vrm30','Tcurv31','Pcurv32','Ele33','Slope34','Asp_c35','Asp_s36','Est37','Nor38','Roug39',
           'TPI40','TRI41','Dx42','Dxx43','Dy44','Dyy45','Maj48','lith.8','lith.24','WTD85','smap','SUMAP','Ref_band1',
           'ref_band2','ref_band3','ref_band4','ref_band5','ref_band6','ref_band7','lon','lat','Sed_thick','Peat',
           'TWI63','ALT64','Permafrost65','Peatclass66','CWDgrow67','ecozone68')
logp <- function(x) ifelse(x == 0, -1, log(x))
d10 <- d1 %>% dplyr::select(socd1m, all_of(covs0), fold_col) %>% rename(SOC = socd1m) %>%
  mutate_all(list(~ replace(., . == c(-9999), NA))) %>%
  filter(SOC != 0 & SOC < 85 & SOC > 0.01 & !is.na(SOC)) %>%
  mutate(NDVI1 = NDVI1 / 10000, EVI2 = EVI2 / 10000, NPP3 = NPP3 / 10000, GPP4 = GPP4 / 1000, PET5 = log(PET5), AI6 = AI6 / 10000, BIO4 = BIO4 / 100,
         BIO12 = logp(BIO12), BIO13 = logp(BIO13), BIO14 = logp(BIO14), BIO15 = logp(BIO15), BIO16 = logp(BIO16), BIO17 = logp(BIO17), BIO18 = logp(BIO18), BIO19 = logp(BIO19),
         Solar26 = Solar26 / 10000, So_seasonal27 = So_seasonal27 / 10000, Vrm30 = Vrm30 * 10^4, Tcurv31 = Tcurv31 * 10^4, Pcurv32 = Pcurv32 * 10^4,
         Ref_band1 = Ref_band1 / 10000, ref_band2 = ref_band2 / 10000, ref_band3 = ref_band3 / 10000, ref_band4 = ref_band4 / 10000, ref_band5 = ref_band5 / 10000, ref_band6 = ref_band6 / 10000, ref_band7 = ref_band7 / 10000,
         Maj48 = factor(Maj48, levels = 1:10), lith.8 = factor(lith.8, levels = 1:7), lith.24 = factor(lith.24, levels = 1:15), Permafrost65 = factor(Permafrost65, levels = 0:4)) %>% na.omit()
d11 <- sf::st_drop_geometry(d10)

## step4 model building

h2o.init(nthreads = 36, max_mem_size = '300G')
h2o.no_progress()
d1_hdf <- as.h2o(d11)
set.seed(1199)
xvars <- setdiff(names(d1_hdf), c('SOC', 'fold_col', 'EVI2', 'lon', 'lat', 'Peat', 'ALT64', 'Peatclass66', 'ecozone68'))
h_auto <- h2o.automl(x = xvars, y = 'SOC', training_frame = d1_hdf, exclude_algos = c('DeepLearning', 'GLM'), max_models = NULL,
                     fold_column = 'fold_col', stopping_metric = 'RMSE', keep_cross_validation_predictions = TRUE, seed = 1807)
lb <- as.data.frame(h2o.get_leaderboard(h_auto, extra_columns = 'ALL'))
saveRDS(lb, file.path(model_dir, 'leaderboard_peat_socd.rds'))
gbm1 <- h2o.get_best_model(h_auto, algorithm = 'GBM', criterion = 'RMSE')
h2o.saveModel(gbm1, model_dir, force = TRUE)

## step5 map prediction

ca_covas$Depth <- 0
data0 <- as.data.frame(as.matrix(ca_covas))
ca_base_new <- data0 %>% mutate_all(list(~ replace(., . == c(-9999), NA))) %>%
  mutate(NDVI1 = NDVI1 / 10000, EVI2 = EVI2 / 10000, NPP3 = NPP3 / 10000, GPP4 = GPP4 / 1000, PET5 = log(PET5), AI6 = AI6 / 10000, BIO4 = BIO4 / 100,
         BIO12 = logp(BIO12), BIO13 = logp(BIO13), BIO14 = logp(BIO14), BIO15 = logp(BIO15), BIO16 = logp(BIO16), BIO17 = logp(BIO17), BIO18 = logp(BIO18), BIO19 = logp(BIO19),
         Solar26 = Solar26 / 10000, So_seasonal27 = So_seasonal27 / 10000, Vrm30 = Vrm30 * 10^4, Tcurv31 = Tcurv31 * 10^4, Pcurv32 = Pcurv32 * 10^4,
         Ref_band1 = Ref_band1 / 10000, ref_band2 = ref_band2 / 10000, ref_band3 = ref_band3 / 10000, ref_band4 = ref_band4 / 10000, ref_band5 = ref_band5 / 10000, ref_band6 = ref_band6 / 10000, ref_band7 = ref_band7 / 10000,
         Maj48 = factor(Maj48, levels = 1:10), lith.8 = factor(lith.8, levels = 1:7), lith.24 = factor(lith.24, levels = 1:15), Permafrost65 = factor(Permafrost65, levels = 0:4),
         lon = lon / 10^6, lat = lat / 10^6, Peat = cut(Peat, breaks = c(0, 5, 30, 60, 100), labels = c(5, 30, 60, 100), include.lowest = TRUE, right = TRUE))
rm(data0); gc()
ind <- rowSums(is.na(ca_base_new)) == 0
pred_df <- ca_base_new[ind, ] %>% mutate(across(where(is.numeric), ~ round(., 4)))
rm(ca_base_new); gc()
test <- as.h2o(pred_df); rm(pred_df); gc()
p <- as.data.frame(h2o.predict(gbm1, test))$predict
p <- pmin(p, 85)
template <- rast(ca_covas[[1]])
v <- rep(NA_real_, ncell(template)); v[which(ind)] <- p
r <- setValues(template, v); names(r) <- 'GBM'
writeRaster(r, file.path(map_dir, 'peat_socd_1m.tif'), overwrite = TRUE)
h2o.shutdown(prompt = FALSE)
