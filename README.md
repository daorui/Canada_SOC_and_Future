# Baseline soil organic carbon map of Canada (0–1 m)

Code and data for the baseline soil organic carbon (SOC) map used in:

> Han, D., et al. Warming drives persistent soil carbon losses in Canada despite increasing plant inputs. *Global Change Biology* (submitted).

The repository contains the workflow used to map upper-metre SOC stocks of Canadian mineral soils and peatlands at 30 arc-second (~1 km) resolution, and the resulting national SOC density map.

## Repository contents

```
code/
  step1_data_cleaning.R                    soil profile compilation, peat/mineral classification, organic-horizon removal, depth harmonisation
  step2_mineral_soc_model.R                mineral SOC concentration model (h2o AutoML) and prediction at 0, 5, 15, 30, 60 and 100 cm
  step3_bulk_density_model.R               mineral bulk density model (h2o AutoML) and prediction at the same depths
  step4_peat_socd_model.R                  peat 0–1 m carbon density model (h2o AutoML, GBM) and prediction
  step5_carbon_stocks_and_national_map.R   mineral and peat carbon density in the upper metre and the national map
data/
  Canada_SOC_0-100cm_kgC_m-2_1km.tif       national SOC density map, 0–1 m
```

Scripts are run in order from the repository root.

## National SOC density map

| Item | Value |
|---|---|
| File | `data/Canada_SOC_0-100cm_kgC_m-2_1km.tif` |
| Variable | Soil organic carbon density, 0–1 m (mineral soil + peat) |
| Unit | kg C m⁻² of grid-cell area |
| Resolution | 30 arc seconds (~1 km) |
| Coordinate reference system | WGS 84 (EPSG:4326) |
| Data type | UInt16 with scale factor 0.01 (values are returned in kg C m⁻² by GDAL-based readers, e.g. `terra::rast()`); precision 0.01 kg C m⁻² |
| NoData | 65535 |
| Domain | Canada, excluding water, snow and ice, barren land and built-up areas (2020 Land Cover of Canada) |
| National stock | 124.6 Pg C (mineral soil 78.8 Pg C; peat 45.7 Pg C) |

Mineral and peat densities are weighted by their area fractions within each grid cell (peat fraction from PEATMAP), so that the two components sum to the total density of the cell.

## Input data

Input data are not redistributed in this repository. They can be obtained from the original sources:

| Input | Source |
|---|---|
| Mineral and peat soil profiles | World Soil Information Service (WoSIS), September 2019 snapshot (Batjes et al., 2020) |
| Peat cores | Canadian peat-profile database (Bauer et al., 2024) |
| Environmental covariates | Listed with sources in Supplementary Table S1 of the paper |
| Peat fraction | PEATMAP (Xu et al., 2018) |
| Peatland classes | Thompson et al. (2016); Webster et al. (2018) |
| Depth to bedrock | Shangguan et al. (2017) |
| Land cover mask | 2020 Land Cover of Canada (Latifovic, 2022) |

The scripts expect the following local files:

```
data/raw/wosis2019_soc_layers.csv
data/raw/wosis2019_bd_layers_canada.csv
data/raw/canadian_peat_profile_database.rds
data/covariates/covariates_mineral.tif
data/covariates/covariates_bd.tif
data/covariates/covariates_peat.tif
data/covariates/peat_fraction.tif
data/covariates/peatland_classes.tif
data/covariates/layer_thickness_to_bedrock_5layers.tif
data/covariates/land_cover_mask.tif
```

## Software

R 4.5.3 with terra 1.9.46, sf 1.1.2, dplyr 1.1.4, tidyr 1.3.1, aqp 1.42, ithir 1.0, spatialsample 0.4.0, blockCV 3.2.0 and h2o 3.40.0.4.

## References

- Batjes, N. H., Ribeiro, E., & van Oostrum, A. (2020). Standardised soil profile data to support global mapping and modelling (WoSIS snapshot 2019). *Earth System Science Data*, 12, 299–320.
- Bauer, I. E., et al. (2024). Peat profile database from peatlands in Canada. *Ecology*, 105, e4398.
- Latifovic, R. (2022). 2020 Land Cover of Canada [Dataset]. Natural Resources Canada.
- Shangguan, W., Hengl, T., Mendes de Jesus, J., Yuan, H., & Dai, Y. (2017). Mapping the global depth to bedrock for land surface modeling. *Journal of Advances in Modeling Earth Systems*, 9, 65–88.
- Thompson, D. K., Simpson, B. N., & Beaudoin, A. (2016). Using forest structure to predict the distribution of treed boreal peatlands in Canada. *Forest Ecology and Management*, 372, 19–27.
- Webster, K. L., et al. (2018). Spatially-integrated estimates of net ecosystem exchange and methane fluxes from Canadian peatlands. *Carbon Balance and Management*, 13, 16.
- Xu, J., Morris, P. J., Liu, J., & Holden, J. (2018). PEATMAP: Refining estimates of global peatland distribution based on a meta-analysis. *CATENA*, 160, 134–140.

## License

- Code (`code/`): MIT License (see `LICENSE`).
- Data (`data/`): Creative Commons Attribution 4.0 International (CC BY 4.0) (see `data/LICENSE.md`).

## Contact

Daorui Han, McMaster University.
