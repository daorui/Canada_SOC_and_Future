library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(spatialsample)
library(h2o)

model_dir <- 'models/mineral_soc/'
map_dir <- 'maps/'
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

## step1 training data and covariate extraction

ca_covas <- rast('data/covariates/covariates_mineral.tif')
soc_mineral <- readRDS('data/processed/mineral_soc_splined.rds')
soc_mineral_sf <- st_as_sf(soc_mineral, coords = c('lon', 'lat'), crs = 4326)
soc_mineral_cov <- terra::extract(ca_covas, soc_mineral_sf)
soc_mineral_u <- cbind(soc_mineral_sf, soc_mineral_cov)

## step2 spatial cross-validation folds

spatialid <- soc_mineral_u %>% group_by(id) %>% reframe(geometry = geometry[1]) %>% st_as_sf()
ca_data_u <- spatialid
set.seed(5234)
folds <- spatial_clustering_cv(ca_data_u, v = 10)
fold_n <- function(x) assessment(x)$id
fold_n1 <- purrr::map(folds$splits, fold_n)
new_group <- 1:nrow(ca_data_u)
for (i in 1:nrow(folds)) new_group[which(ca_data_u$id %in% fold_n1[[i]])] <- i
ca_data_u$fold_col <- new_group
soc_mineral_u2 <- st_join(soc_mineral_u, ca_data_u[, c(1, 3)])

## step3 covariate preprocessing

d1 <- soc_mineral_u2
item0 <- paste0('X', c(5, 15, 30, 45, 60, 80, 100), 'cm')
names(d1)[2:8] <- item0
covs0 <- c('NDVI1','EVI2','NPP3','GPP4','PET5','AI6','BIO1','BIO10','BIO11','BIO12','BIO13','BIO14','BIO15','BIO16',
           'BIO17','BIO18','BIO19','BIO2','BIO3','BIO4','BIO5','BIO6','BIO7','BIO8','BIO9','Solar26','So_seasonal27',
           'Pop_den28','HFP29','Vrm30','Tcurv31','Pcurv32','Ele33','Slope34','Asp_c35','Asp_s36','Est37','Nor38','Roug39',
           'TPI40','TRI41','Dx42','Dxx43','Dy44','Dyy45','Maj48','lith.8','lith.24','WTD85','smap','SUMAP','Ref_band1',
           'ref_band2','ref_band3','ref_band4','ref_band5','ref_band6','ref_band7','lon','lat','Sed_thick','Peat',
           'TWI63','ALT64','Permafrost65','Peatclass66','CWD67','ecozone68')
logp <- function(x) ifelse(x == 0, -1, log(x))
d10 <- d1 %>% dplyr::select(all_of(item0), all_of(covs0), fold_col) %>%
  tidyr::pivot_longer(cols = all_of(item0), names_to = 'Depth', names_transform = list(Depth = function(x) as.integer(gsub('X|cm', '', x))),
                      values_drop_na = TRUE, values_to = 'SOC') %>%
  mutate(Depth = case_when(Depth == 5 ~ 2.5, Depth == 15 ~ 12.5, Depth == 30 ~ 22.5, Depth == 45 ~ 37.5, Depth == 60 ~ 52.5, Depth == 80 ~ 70, Depth == 100 ~ 90)) %>%
  mutate_all(list(~ replace(., . == c(-9999), NA))) %>%
  filter(SOC != 0 & SOC < 580 & SOC > 0.01 & !is.na(SOC)) %>%
  mutate(SOC = log(SOC * 100)^2,
         NDVI1 = NDVI1 / 10000, EVI2 = EVI2 / 10000, NPP3 = NPP3 / 10000, GPP4 = GPP4 / 1000, PET5 = log(PET5), AI6 = AI6 / 10000, BIO4 = BIO4 / 100,
         BIO12 = logp(BIO12), BIO13 = logp(BIO13), BIO14 = logp(BIO14), BIO15 = logp(BIO15), BIO16 = logp(BIO16), BIO17 = logp(BIO17), BIO18 = logp(BIO18), BIO19 = logp(BIO19),
         Solar26 = Solar26 / 10000, So_seasonal27 = So_seasonal27 / 10000, Vrm30 = Vrm30 * 10^4, Tcurv31 = Tcurv31 * 10^4, Pcurv32 = Pcurv32 * 10^4,
         Ref_band1 = Ref_band1 / 10000, ref_band2 = ref_band2 / 10000, ref_band3 = ref_band3 / 10000, ref_band4 = ref_band4 / 10000, ref_band5 = ref_band5 / 10000, ref_band6 = ref_band6 / 10000, ref_band7 = ref_band7 / 10000,
         Maj48 = factor(Maj48, levels = 1:10), lith.8 = factor(lith.8, levels = 1:7), lith.24 = factor(lith.24, levels = 1:15),
         ecozone68 = factor(ecozone68, levels = c(1:12, 14:31))) %>% na.omit()
d11 <- sf::st_drop_geometry(d10) %>% filter(Depth < 110)

## step4 model building

h2o.init(nthreads = 36, max_mem_size = '300G')
h2o.no_progress()
d1_hdf <- as.h2o(d11)
set.seed(1199)
xvars <- setdiff(names(d1_hdf), c('SOC', 'fold_col', 'NDVI1', 'EVI2', 'lon', 'lat', 'Peat', 'ALT64', 'Peatclass66', 'Permafrost65', 'CWD67'))
h_auto <- h2o.automl(x = xvars, y = 'SOC', training_frame = d1_hdf, exclude_algos = c('DeepLearning', 'GLM'), max_models = NULL,
                     fold_column = 'fold_col', stopping_metric = 'RMSE', keep_cross_validation_predictions = TRUE, seed = 1807)
lb <- as.data.frame(h2o.get_leaderboard(h_auto, extra_columns = 'ALL'))
saveRDS(lb, file.path(model_dir, 'leaderboard_mineral_soc.rds'))
stack_m <- h2o.get_best_model(h_auto, algorithm = 'StackedEnsemble', criterion = 'RMSE')
h2o.saveModel(stack_m, model_dir, force = TRUE)

## step5 map prediction

ca_covas$Depth <- 0
data0 <- as.data.frame(as.matrix(ca_covas))
ca_base_new <- data0 %>% mutate_all(list(~ replace(., . == c(-9999), NA))) %>%
  mutate(NDVI1 = NDVI1 / 10000, EVI2 = EVI2 / 10000, NPP3 = NPP3 / 10000, GPP4 = GPP4 / 1000, PET5 = log(PET5), AI6 = AI6 / 10000, BIO4 = BIO4 / 100,
         BIO12 = logp(BIO12), BIO13 = logp(BIO13), BIO14 = logp(BIO14), BIO15 = logp(BIO15), BIO16 = logp(BIO16), BIO17 = logp(BIO17), BIO18 = logp(BIO18), BIO19 = logp(BIO19),
         Solar26 = Solar26 / 10000, So_seasonal27 = So_seasonal27 / 10000, Vrm30 = Vrm30 * 10^4, Tcurv31 = Tcurv31 * 10^4, Pcurv32 = Pcurv32 * 10^4,
         Ref_band1 = Ref_band1 / 10000, ref_band2 = ref_band2 / 10000, ref_band3 = ref_band3 / 10000, ref_band4 = ref_band4 / 10000, ref_band5 = ref_band5 / 10000, ref_band6 = ref_band6 / 10000, ref_band7 = ref_band7 / 10000,
         Maj48 = factor(Maj48, levels = 1:10), lith.8 = factor(lith.8, levels = 1:7), lith.24 = factor(lith.24, levels = 1:15), ecozone68 = factor(ecozone68, levels = c(1:12, 14:31)),
         lon = lon / 10^6, lat = lat / 10^6, Peat = cut(Peat, breaks = c(0, 5, 30, 60, 100), labels = c(5, 30, 60, 100), include.lowest = TRUE, right = TRUE))
rm(data0); gc()
ind <- rowSums(is.na(ca_base_new)) == 0
pred_df <- ca_base_new[ind, ] %>% mutate(across(where(is.numeric), ~ round(., 4)))
rm(ca_base_new); gc()
test <- as.h2o(pred_df); rm(pred_df); gc()
template <- rast(ca_covas[[1]])
for (i in c(0, 5, 15, 30, 60, 100)) {
  test[, 'Depth'] <- i
  p <- as.data.frame(h2o.predict(stack_m, test))$predict
  p <- pmin(p, 120.3014); p <- exp(sqrt(p)) / 100
  v <- rep(NA_real_, ncell(template)); v[which(ind)] <- p
  r <- setValues(template, v); names(r) <- 'Stack'
  writeRaster(r, file.path(map_dir, paste0('mineral_soc_', i, 'cm.tif')), overwrite = TRUE)
}
h2o.shutdown(prompt = FALSE)
