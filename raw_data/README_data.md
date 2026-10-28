# Data description

## Input data (`raw_data/`)

### `Baskiera_Gvozdik_Dataset_B.txt`

Tab-delimited. Annual measurements of 37 alpine newts (*Ichthyosaura alpestris*; 18 females, 19 males) from a population near Jihlava, Czech Republic (49.2247 N, 16.0579 E; 600 m a.s.l.), monitored from metamorphosis (2014) to 2018. One row per individual and measurement (178 rows). Originally published in Baškiera & Gvoždík (2021), *Comp. Biochem. Physiol. A*, 253, 110858.

| Column | Description | Units |
|---|---|---|
| `Year` | Year of measurement | year |
| `date` | Date of measurement (DD.MM.YYYY) | date |
| `ID` | Individual identifier (matches the folder names in `ind/`) | — |
| `Sex` | Sex: `F` female, `M` male | — |
| `BM` | Body mass | g |
| `VO2` | Standard (resting) metabolic rate, as oxygen consumption | mL O2 h⁻¹ |
| `resVO2` | Residual of VO2 after accounting for body mass, as in the original study (not used in the analyses) <!-- CONFIRM --> | mL O2 h⁻¹ |

### `Females2014_clutchsize_v2.csv`

Comma-separated. Clutch data for 19 females of the same cohort (2018 breeding season).

| Column | Description | Units |
|---|---|---|
| `ID` | Individual identifier | — |
| `Number of eggs` | Number of eggs laid in 2018 | eggs |
| `Body mass (6.4.2018)` | Body mass on 6 April 2018 | g |
| `Start_laying` | Date the female started laying in 2018 (M/D/YYYY) | date |
| `First laying` | Year of first reproduction | year |

### `alpestris_shapefile/`

**Not included** (IUCN Red List terms of use do not allow redistribution). Polygon of the current distribution of *I. alpestris* (IUCN 2021, https://doi.org/10.2305/IUCN.UK.2022-1.RLTS.T59472A89702213.en), used only to outline the range on the maps. See the main `README.md` for download instructions.

## Derived outputs (`results/`, archived as `results.zip`)

All metric columns use the following names and units (lifetime totals are summed over each simulated life):

| Metric | Description | Units |
|---|---|---|
| `MaxWgt` | Maximum body mass | g (wet) |
| `MaxLen` | Maximum snout–vent length | cm |
| `Tmax`, `Tmin` | Maximum and minimum annual body temperature | °C |
| `MinRes` | Minimum reserve energy | J |
| `Bsk` | Lifetime basking time (above ground without foraging) | h |
| `Forage` | Lifetime foraging time | h |
| `Food` | Lifetime food intake | g (dry) |
| `Fecundity` | Lifetime number of eggs (females) | eggs |
| `ReproEnergy` | Lifetime reproductive energy, `Fecundity × E_0` (used for males) | J |
| `Clutches` | Lifetime number of clutches | count |
| `LifeSpan` | Lifespan | years |
| `R0` | Net reproductive rate | — |
| `MeanDailyDef` | Mean daily water need (`TotDeficit / DaysDeficit`) | g day⁻¹ |
| `TotDeficit` | Total water deficit over days with negative water balance | g |
| `DaysDeficit` | Number of days with negative water balance | days |
| `Pct_Des` | Maximum desiccation level reached (NicheMapR output) | % |
| `CauseDeath` | NicheMapR cause-of-death code (0 = survived the whole simulation; see NicheMapR `ectotherm()` documentation) | code |

| File | Rows | Content |
|---|---|---|
| `europe_results{,_w2,_w4}.rds` | one per grid cell × parameterisation | Raw simulation output per 10′ grid cell for the current climate, +2 °C (`_w2`) and +4 °C (`_w4`). Columns: `type` (`species` = population-average parameterisation, `individual`), `id` (individual ID, `sp` for population average), `sex`, `lon`, `lat` (decimal degrees, WGS84), `elev` (m), and the metrics above. |
| `europe_cleaned{,_w2,_w4}.rds` | as above | As `europe_results`, with error rows removed and added columns `survived`, `dehydrated`, `temp_stress` (logical flags derived from `CauseDeath`; see `R/4.1_analyze_europe.R`), `E_0` (initial egg energy, J) and `ReproEnergy`. |
| `europe_summary{,_w2,_w4}.rds` | one per grid cell × sex | Per-cell summary across individuals: `<metric>_mean`, `_median`, `_sd`, `_p5`, `_p95`, `_min`, `_max`, `_cv`; `frac_survived`, `frac_dehydrated`, `frac_temp_stress` (fraction of individuals); `sp_<metric>` (population-average prediction). |
| `europe_jensen{,_w2,_w4}.rds` | one per grid cell × sex × metric | `P_species` (population-average prediction, P_pop), `P_ind` (mean individual prediction), `J = P_species − P_ind` (metric units), `J_rel` (relative bias), `sign` (sign of J). |
| `jensen_comparison.rds` | list | `per_cell`: J, |J|, predictions and sign for current (`_cur`), +2 °C (`_w2`) and +4 °C (`_w4`), with `delta_abs_J_w2/_w4` = |J_future| − |J_current|, sign changes and direction of change. `per_metric`: per sex and metric summaries (mean |J|, mean/median Δ|J|, % cells with increase/decrease/sign change, relative Δ|J| in %). |
| `param_cv.rds` | 4 values | Coefficient of variation (%) of the individually estimated DEB parameters (`p_M`, `v`, `z_m`, `z`). |
| `mat_europe.tif` | raster | Mean annual air temperature (°C) from New et al. (2002), cropped to Europe. |
