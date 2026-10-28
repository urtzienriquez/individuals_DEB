# Beyond the average phenotype: individual variation alters mechanistic forecasts of climate responses

Data and code for:

> Enriquez-Urzelai, U., Pascucci, G. & Gvoždík, L. Beyond the average phenotype: individual variation alters mechanistic forecasts of climate responses. Submitted to *Ecology Letters*.

Archived version (data + code): <https://doi.org/10.5281/zenodo.XXXXXXX>

The study estimates Dynamic Energy Budget (DEB) parameters for a population-average alpine newt (*Ichthyosaura alpestris*) and for 37 individually monitored newts, couples each parameterisation with NicheMapR biophysical simulations across Europe (current climate, +2 °C and +4 °C), and quantifies Jensen bias, i.e. the difference between the prediction of the average phenotype and the mean of the individual predictions.

## Repository layout

| Path | Content |
|---|---|
| `raw_data/` | Input data: individual body mass and metabolic rate (`Baskiera_Gvozdik_Dataset_B.txt`) and clutch sizes (`Females2014_clutchsize_v2.csv`). Column descriptions in [`raw_data/README_data.md`](raw_data/README_data.md). The species range polygon (`alpestris_shapefile/`) is **not included** (IUCN terms of use); see below. |
| `sp/` | AmP/DEBtool (MATLAB) files for the population-average DEB parameter estimation (`mydata_`, `pars_init_`, `predict_`, `run_` files) and its output (`results_*.mat`, html reports). `sp_keep/` is a backup copy of the final estimation. |
| `ind/<ID>/` | The same AmP/DEBtool files for each individual (folder name = individual ID in the raw data). |
| `R/` | R scripts (run order below) and shared functions (`utils.R`, `params_definition.R`, `parse_ind_pars.R`). |
| `results/` | Derived simulation outputs (`.rds`, `mat_europe.tif`). **Not in the git repository because of their size; included in the Zenodo archive as `results.zip`.** Described in `raw_data/README_data.md`. |
| `figures/1_base/` | Base figures (PDF) written by the R scripts. `figures/2_final/` holds standalone LaTeX/TikZ files that add labels to the base figures; compile each with `pdflatex` to get the final figures. |
| `manuscript/` | Manuscript sources (knitr `.Rnw` + LaTeX). `1_index.Rnw` builds the main text and Appendices S1–S3; all numbers in the text are computed from `results/` by `R/source_results_manuscript.R`. |

## Software

- R 4.5.0
- NicheMapR **3.3.3** (https://github.com/mrke/NicheMapR; the exact version matters for reproducing the simulations)
- R packages: terra 1.9.11, sf 1.1.1, mgcv 1.9.4, dplyr 1.2.1, tidyr, purrr, ggplot2 4.0.3, maps, mapdata, tikzDevice, R.matlab, xml2, future, future.apply
- MATLAB with DEBtool_M and AmPtool (https://github.com/add-my-pet) for the DEB parameter estimation
- LaTeX (pdflatex, biber) and the `knitrmini` R package to build the manuscript and figure labels

## Reproducing the analyses

All R scripts are run from the `R/` directory (paths are relative, e.g. `../results/`).

### Shortcut: reproduce figures and numbers from archived results (no simulations)

1. Download `results.zip` from Zenodo and unzip it into `results/`. Download the range polygon (see *Range polygon* below) into `raw_data/alpestris_shapefile/`.
2. Run the figure scripts: `4.2_map_europe.R` (Fig. 2, Figs S3.1–S3.3), `7.1_map_plus_curve.R` (Fig. 3, Fig. S3.4) and `7.2_map_suitability.R` (Fig. 4, Figs S3.5–S3.7). They write base PDFs to `figures/1_base/`.
3. Compile the label files in `figures/2_final/` with `pdflatex`.
4. Knit the manuscript: in `manuscript/`, run `Rscript -e 'knitrmini::knit2pdf("1_index.Rnw")'`. Every number reported in the text is recomputed from `results/` at this step.

### Full pipeline

| Step | Script(s) | Needs | Output |
|---|---|---|---|
| 1. DEB estimation, population average | `sp/run_Ichthyosaura_alpestris.m` | MATLAB + DEBtool/AmPtool | `sp/results_*.mat`, html report |
| 1. DEB estimation, individuals | `1.1_setup_ind_templates.R`, `1.2_prepare_ind_files.R` (write the `ind/<ID>/` files from the raw data), then `ind/<ID>/run_*.m` | R; MATLAB | `ind/<ID>/results_*.mat` |
| 2. Parameter analysis | `2.1_explore_pars.R`, `2.2_analyze_pars.R`, `2.3_growth_plots.R` | R (2.3 runs NicheMapR) | `results/param_cv.rds`, Figs S2.x |
| 3. European simulations | `3.0_submit_{baseline,w2,w4}.sh` → `3.1_i_europe_chunk.sh` → `3.1_run_europe_chunk.R`, then `3.2_merge_chunks.R` | HPC (PBS), NicheMapR | `results/europe_results{,_w2,_w4}.rds` |
| 4. Current-climate post-processing | `4.1_analyze_europe.R`, `4.2_map_europe.R` | R | `europe_{cleaned,summary,jensen}.rds`, figures |
| 5. Warming post-processing | `5.1_analyze_warming.R` (`5.2_map_warming.R` maps not used in the paper) | R | `europe_*_{w2,w4}.rds` |
| 6. Comparison across scenarios | `6.1_compare_current_future.R` (`6.2`, `6.3` not used in the paper) | R | `results/jensen_comparison.rds` |
| 7. Final figures | `7.1_map_plus_curve.R`, `7.2_map_suitability.R` | R | Figs 3, 4, S3.4–S3.7 |

Notes:
- The European simulations take roughly 30 chunks × 12 h on 10 cores per climate scenario (see `3.1_i_europe_chunk.sh`).
- `7.1_map_plus_curve.R` builds `results/mat_europe.tif` (mean annual air temperature) from the New et al. (2002) climatology distributed with NicheMapR (`global_climate.nc`, obtained with `NicheMapR::get.global.climate()`) if the file is missing; the archived file is used otherwise.
- `example_species.R` and `example_individual_8045.R` run single-site simulations for illustration.

## Data sources

- Individual growth and metabolic data: Baškiera, S. & Gvoždík, L. (2021). Repeatability and heritability of resting metabolic rate in a long-lived amphibian. *Comp. Biochem. Physiol. A Mol. Integr. Physiol.*, 253, 110858. https://doi.org/10.1016/j.cbpa.2020.110858
- Species-level life-history data: Gvoždík & Kristín (2017); Miaud et al. (2000) (see manuscript).
- Climate: New, M., Lister, D., Hulme, M. & Makin, I. (2002). A high-resolution data set of surface climate over global land areas. *Clim. Res.*, 21, 1–25.
- Range polygon: IUCN (2021), *Ichthyosaura alpestris*, The IUCN Red List of Threatened Species, https://doi.org/10.2305/IUCN.UK.2022-1.RLTS.T59472A89702213.en

### Range polygon

The IUCN range map cannot be redistributed. It is used only to draw the range outline on the maps (scripts `4.2`, `5.2`, `6.2`, `6.3`, `7.1`, `7.2`) and does not affect any analysis. To reproduce the outlines, download the *Ichthyosaura alpestris* range polygon from the IUCN Red List spatial data portal (https://www.iucnredlist.org/resources/spatial-data-download; free registration required), and save it as `raw_data/alpestris_shapefile/data_0.shp` (with its `.dbf`, `.shx`, `.prj` and `.cpg` companion files).

## License

Code: MIT (see `LICENSE`). Data and derived outputs: CC BY 4.0.
