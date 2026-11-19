library(terra)

map_dir <- 'maps/'
cov_dir <- 'data/covariates/'
out_dir <- 'output/'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## step1 mineral soil carbon density in the upper metre

deps <- c(0, 5, 15, 30, 60, 100)
soc_d <- list()
bd_d <- list()
for (k in 1:6) {
  soc_d[[k]] <- rast(file.path(map_dir, paste0('mineral_soc_', deps[k], 'cm.tif')))
  bd_d[[k]] <- rast(file.path(map_dir, paste0('mineral_bd_', deps[k], 'cm.tif')))
}
soc_d <- rast(soc_d)
bd_d <- rast(bd_d)
all_maps4 <- (soc_d[[1:5]] + soc_d[[2:6]]) / 2
f5 <- (bd_d[[1:5]] + bd_d[[2:6]]) / 2
depth00 <- rast(file.path(cov_dir, 'layer_thickness_to_bedrock_5layers.tif'))
peat_frac <- rast(file.path(cov_dir, 'peat_fraction.tif'))
ca_mineral_top <- sum(all_maps4[[1:3]] * f5[[1:3]] * depth00[[1:3]]) * (1 - peat_frac / 100)
ca_mineral_sub <- sum(all_maps4[[4:5]] * f5[[4:5]] * depth00[[4:5]]) * (1 - peat_frac / 100)
top <- ca_mineral_top / 100
sub <- ca_mineral_sub / 100

## step2 peat carbon density in the upper metre

ca_peat <- rast(file.path(map_dir, 'peat_socd_1m.tif'))
peat_type <- rast(file.path(cov_dir, 'peatland_classes.tif'))
peat_ratio_u <- data.frame(peatclass = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                           ratio = c(0.276, 0.2986, 0.2986, 0.2986, 0.2557, 0.2149, 0.27, 0.2638, 0.2667, 0.27))
peat_map_ratio <- classify(peat_type, rcl = as.matrix(peat_ratio_u))
peat_socd_acrotelm <- ca_peat * (peat_frac / 100) * peat_map_ratio
peat_socd_catotelm <- ca_peat * (peat_frac / 100) - peat_socd_acrotelm

## step3 national soil organic carbon density map

mn <- top + sub
pt <- peat_socd_acrotelm + peat_socd_catotelm
lc <- rast(file.path(cov_dir, 'land_cover_mask.tif'))
tot <- ifel(is.na(mn) & is.na(pt), NA, ifel(is.na(mn), 0, mn) + ifel(is.na(pt), 0, pt))
tot <- mask(tot, lc)
names(tot) <- 'SOC_0_100cm_kgCm2'
writeRaster(tot, file.path(out_dir, 'Canada_SOC_0-100cm_kgC_m-2.tif'), overwrite = TRUE,
            gdal = c('COMPRESS=DEFLATE', 'PREDICTOR=3', 'ZLEVEL=6'), datatype = 'FLT4S')
