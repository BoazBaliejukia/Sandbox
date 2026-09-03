# An annually updated event-niche pipeline for Ebolavirus spillover in equatorial Africa: methods from KSPH Code and a documented 2001–2026 presence linelist

**Manuscript type.** Methods and data note (unpublished draft).  
**Provenance.** This draft documents the public “Event Ecological Niche Pipeline” in [carsontelford/KSPH-Code](https://github.com/carsontelford/KSPH-Code) and the 3 September 2026 snapshot of the fork [BoazBaliejukia/KSPH_Code](https://github.com/BoazBaliejukia/KSPH_Code). The pipeline is a modular successor to the presence–background boosted-regression-tree (BRT) workflow published by Telford et al. (2025) [1]. Fitted ensembles, Earth Engine exports, and prediction rasters are gitignored and were **not** present in the public repository at the time of writing. No model-performance numbers, relative-odds maps, or TreeSHAP values are reported as if they had been computed in this run.

---

## Abstract

**Background.** Zoonotic spillover of Ebolavirus species into humans remains rare but operationally important. Forest cover, forest loss, fragmentation, human population, and climate have been used to rank equatorial African landscapes by relative odds of an index case. Telford et al. (2025) fitted annually updated BRT ensembles for 2001–2021. Three subsequent outbreaks—Sudan virus disease in Uganda (2025), *Orthoebolavirus zairense* in Bulape, Democratic Republic of the Congo (2025), and Bundibugyo virus disease centred on Ituri Province (2026)—extend the period that a refreshable pipeline must cover.

**Objective.** To describe the KSPH Code event-niche pipeline from the committed scripts and configuration files, and to tabulate the deposited presence linelist without fabricating covariate extracts, model metrics, or maps that the public repository does not contain.

**Methods.** We reviewed `R_python_code/` and `config/` in the BoazBaliejukia fork (identical in structure to the upstream Telford repository). Presence rows were read from `config/outcomes_csv.csv`. Year, species-code, country, and geographic ranges were tabulated directly from that file. Script constants were recorded for pseudo-absence sampling, Google Earth Engine products, BRT ensembles, temporal forward validation, SuperLearner, and relative-odds-ratio (ROR) construction. We did not authenticate Earth Engine, fit models, or create `data/`, `models/`, or `outputs/` artefacts.

**Results.** The linelist contains 27 presence records dated 2001–2026: 18 coded Z (Zaire), 6 coded S (Sudan), and 3 coded B (Bundibugyo), spanning Gabon, the Republic of the Congo, South Sudan, the Democratic Republic of the Congo (DRC), Uganda, and Guinea. Coordinates range from 5.26°S to 8.62°N and from 10.06°W to 32.73°E. Script `01_make_dataset1.R` retains years 2001–2025 only, so the 2026 Bunia row is excluded from the intended training table; the remaining 26 points all fall inside the configured equatorial bounding box. The pipeline is designed to sample 5,000 random land pseudo-absences, extract 35 named predictors over donut buffers via Hansen v1.13, LandScan, ERA5-Land, MODIS NDVI, and SRTM, fit a 50-member Bernoulli BRT ensemble (`dismo::gbm.step`; tree complexity 5; learning rate 0.0012) with 50 controls per event, and write 2020–2025 ROR rasters plus a 2022–2025 temporal-forward-validation series. Those extracts, models, and rasters are not in the public repository.

**Conclusions.** KSPH Code is a handoff-ready specification for annually updated event-niche modelling, and the committed linelist documents three outbreaks after the Telford et al. training window. It is not, by itself, a completed analysis: independent Results on discrimination, variable importance, or where 2025–2026 events sat on predicted ROR surfaces still require a local Earth Engine run.

**Keywords.** Ebolavirus; zoonotic spillover; ecological niche model; boosted regression trees; SuperLearner; equatorial Africa; open methods.

---

## Outline (IMRaD)

| Section | Content |
|---|---|
| **Title** | An annually updated event-niche pipeline for Ebolavirus spillover in equatorial Africa |
| **Abstract** | Background, aim, script-based methods, linelist-only results, limitation that models were not run |
| **Introduction** | Spillover ecology; gap between static niche maps and annually refreshable odds; 2025–2026 outbreaks; rationale for documenting this pipeline |
| **Methods** | Study design; linelist; study area; GEE covariates; BRT, temporal validation, SuperLearner; ROR definition; what this draft did and did not run |
| **Results** | Composition of `outcomes_csv.csv`; year-filter implication for 2026; intended products that are missing from Git |
| **Discussion** | Comparison with Telford et al. and earlier niche studies; implications of 2025–2026 events; limitations of presence–background designs and of an unrun pipeline |
| **Conclusion** | The repository is a reproducible recipe plus an updated linelist, not a finished predictive product |

---

## 1. Introduction

Outbreaks of Ebola virus disease typically begin when an Ebolavirus species infects a person from an animal source and then spreads in human networks [1,2]. Four species have caused human disease: *Zaire ebolavirus* (*Orthoebolavirus zairense*, EBOV), *Sudan ebolavirus* (*Orthoebolavirus sudanense*, SUDV), *Bundibugyo ebolavirus* (*Orthoebolavirus bundibugyoense*, BDBV), and *Taï Forest ebolavirus* [1]. A definitive reservoir has not been established for each species, although forest-associated bats are leading candidates and field investigations have implicated fruit bats in at least some transmission settings [3–5]. Index-case investigations have also implicated hunting, butchering, and market handling of wildlife [2,5].

The public-health burden of secondary transmission can be extreme. The 2014–2016 West African epidemic caused nearly 30,000 reported cases and more than 11,000 deaths, and the 2018–2020 outbreak in the DRC caused approximately 3,500 cases and 2,300 deaths [1]. Identification of spillover events has become more frequent in the twenty-first century than in the first two decades after 1976, reflecting some combination of improved detection and changing human–animal contact [1,6]. Because spillover remains rare in any given place and year, surveillance cannot be equally intensive everywhere. Models that rank locations by relative odds of an index case, using covariates that are updated annually, are therefore of operational interest [1,7].

Several ecological studies have associated EVD emergence with forest structure and climate. Pigott and colleagues mapped a zoonotic niche for Ebola virus disease using environmental suitability [8]. Schmidt et al. related human and great-ape spillover intensity to vegetation, population density, and rainfall seasonality [9]. Rulli et al. found that index cases during 2004–2014 occurred preferentially in fragmentation hotspots rather than in hotspots of forest *loss* per se [10]. Olivero et al. reported that outbreaks along rainforest margins were associated with closed-forest loss in the preceding two years [11]. Those results motivate models that treat forest loss and fragmentation as time-varying predictors, not only as static habitat layers.

Telford et al. addressed a practical gap: high-resolution forest-change and population datasets are released annually, typically with a lag of several months, and can be used to refresh relative-odds surfaces rather than to produce a single historical suitability map [1,12,13]. They fitted BRT ensembles for 2001–2021 and generated relative odds ratios for 2021 and 2022 [1]. Code and a cleaned presence–absence table were deposited separately [14]. BRT methods combine regression trees with boosting and are widely used for ecological prediction when many correlated predictors and nonlinear effects are expected [15]. They are not a substitute for etiologic identification of a reservoir or a causal decomposition of land-use effects [1,16].

Three outbreaks after that paper’s 2022 hold-out year change the operational question. In January 2025, Uganda declared an outbreak of Sudan virus disease; 14 cases (12 confirmed, 2 probable) and four deaths were reported from seven districts, including Kampala and Wakiso, before the outbreak was declared over on 26 April 2025 [17]. In September 2025, the DRC declared an outbreak of Ebola virus disease in Bulape Health Zone, Kasai Province—the country’s 16th EVD outbreak—which reached 64 cases (53 confirmed, 11 probable) and 45 deaths (case-fatality ratio 70.3%) and was declared over on 1 December 2025 [18,19]. In May 2026, laboratory confirmation of Bundibugyo virus in Ituri Province opened the DRC’s 17th Ebola-disease outbreak; WHO reported that the event was thought to have originated in Mongbwalu Health Zone, with subsequent cases in Rwampara and Bunia, and later described large-scale expansion beyond Ituri [20,21]. Those events create a concrete need for a pipeline that can ingest an updated linelist, rebuild annual covariates through at least 2025, and hold recent years out of training.

The KSPH Code repository is that pipeline. It is written as a general “event ecological niche” workflow—presence sampling, Earth Engine extraction, BRT prediction, temporal forward validation, SuperLearner, and interactive reports—while the committed outcome file is an Ebolavirus linelist that now includes 2025 and 2026 rows. The present draft has two aims: (i) to write Methods from the actual scripts and configuration files, and (ii) to report only those Results that can be computed from files that Git tracks, while stating explicitly what the code is designed to produce but has not deposited.

---

## 2. Methods

### 2.1 Study design and intended outcome

The analysis specified in the repository is a presence–background (presence–pseudo-absence) predictive study of documented Ebolavirus events. Script `01_make_dataset1.R` codes `outcome = 1` for rows read from the presence CSV and `outcome = 0` for randomly sampled background points. Script `03_train_predict_brt_simple.R` treats that binary column as a Bernoulli response. An optional `TRAINING_TYPE_FILTER` can restrict events to a single type code (for example `"Z"`) while retaining all control rows; the committed default is an all-types model.

The intended map product is not a probability calibrated to absolute annual incidence. After rasterising ensemble-mean predicted probabilities, the script converts cell probabilities \(p\) to odds \(p/(1-p)\) (with \(p\) bounded away from 0 and 1 by \(10^{-6}\)) and divides by the study-area mean odds in that year, producing a relative odds ratio (ROR) surface, following the definition used by Telford et al. [1]. Derived layers include the 99th-percentile (top 1%) of ROR, one-year ROR ratios, one-year raw-odds ratios, and pixel-relative prediction ratios.

### 2.2 Study area

The README and both R and Python configuration blocks describe mainland Africa within approximately 10° of the Equator. The canonical land mask is `config/africacountries_nolakes.shp`. Script `01_make_dataset1.R` clips that polygon to the bounding box \(\lambda \in [-15.5^\circ, 51.0^\circ]\), \(\varphi \in [-10.0^\circ, 10.0^\circ]\) (`STUDY_AREA_BBOX`). The same box is the Python default (`DEFAULT_STUDY_AREA_BBOX`). A 1° latitude margin is added for covariate extraction so that 100 km-scale donut buffers near the study-area edge can still summarise nearby land. `ALLOW_BBOX_FALLBACK` is `FALSE` by default because sampling the rectangle without the land shapefile can place points in the ocean.

### 2.3 Presence linelist

Presence records are committed at `config/outcomes_csv.csv` (the README still mentions a local `data/outcomes_csv.csv`; the script path used here is the config copy). Required fields after standardisation are longitude, latitude, and year; optional fields are `id`, `type`, and `country`. Type codes in the file are `Z`, `S`, and `B`. The script keeps rows whose year is in `ABSENCE_YEARS` (`2001:2025`) and, when `sf` is available, drops points outside the active study polygon (`FILTER_PRESENCES_TO_STUDY_AREA <- TRUE`).

We treated the CSV as the analysis-ready linelist for this draft. Year, type, country, and coordinate ranges were computed from the 27 data rows. We did not independently re-geocode village names. For the three post-2022 rows we compared identifiers and coordinates with WHO Disease Outbreak News only to describe what the points represent, not to edit the file [17,18,20].

### 2.4 Pseudo-absence sampling

With `sf` installed, `01_make_dataset1.R` samples `N_PSEUDO_ABSENCE = 5000` points uniformly in the bounding box of the cropped land polygon and retains only those that intersect the polygon. Years are drawn with replacement from 2001–2025. The random seed is `20260813`; the sampling version string is `pseudo_absence_random_v1`. `EXCLUDE_AROUND_PRESENCES_M` is 0, so background points are not buffered away from known events. Output is `data/dataset1.csv` (gitignored). Unlike Telford et al., who used LandScan-weighted absences [1,13], this first step is spatially random on land.

### 2.5 Covariate extraction (Google Earth Engine)

Notebook `02_predGrid_trainSet_extraction.ipynb` and helpers `02_drivehelper_covariates.py` and `02_prediction_grid_helper.py` orchestrate Earth Engine exports to Google Drive (`Event_ENM_Exports`). Python dependencies are `earthengine-api`, `pandas`, `geopandas`, and `shapely`.

**Products named in the helper.**

| Product | Earth Engine ID | Role |
|---|---|---|
| Hansen Global Forest Change v1.13 (2000–2025) | `UMD/hansen/global_forest_change_2025_v1_13` | Annual forest cover, same-year and lagged loss, edge [12] |
| LandScan Global | `projects/sat-io/open-datasets/ORNL/LANDSCAN_GLOBAL` | Population density; latest available year fixed at 2023 [13] |
| ERA5-Land monthly aggregates | `ECMWF/ERA5_LAND/MONTHLY_AGGR` | Annual precipitation and 2 m temperature [22] |
| ERA5-Land daily aggregates | `ECMWF/ERA5_LAND/DAILY_AGGR` | Annual potential evapotranspiration (sign-flipped ×1000 to millimetres) [22] |
| MODIS NDVI | `MODIS/061/MOD13A2` | Annual NDVI |
| SRTM elevation | `USGS/SRTMGL1_003` | Elevation, extracted at 5 km |

Hansen-derived variables are coarsened before buffer extraction because native 30 m pixels are too expensive for repeated 100 km summaries (`HANSEN_BUFFER_SCALE_M = 1000`; `HANSEN_REDUCE_MAX_PIXELS = 65535`). Forest cover in a given year is 2000 tree cover with subsequent loss pixels set to 0. Same-year, one-year-prior, and two-year-prior loss indicators are built from the Hansen `lossyear` band; one-year lag is undefined before 2002 and two-year lag before 2003. Fragmentation is the proportion of coarse pixels that are dense forest (cover ≥ 0.40) and have a non-dense 3×3 neighbour (`focal_min`). Climate and NDVI layers include annual totals or means plus anomaly and standardised (z) anomaly relative to a 2000–2024 baseline.

**Buffers.** Training extraction uses four approximate-degree donut rings labelled `0_10km`, `10_25km`, `25_50km`, and `50_100km`, with outer radii 0.090°, 0.225°, 0.449°, and 0.674° (the last label is retained while the outer radius is about 75 km at the Equator). Scaled families (forest cover, forest loss, fragmentation, population) are extracted on all rings; climate, NDVI, and elevation are extracted on `0_10km` only. The prediction-grid helper currently implements the three inner rings. The 35 `include_in_model = true` rows in `config/predictor_list.csv` are the intended BRT/SuperLearner predictors. `config/external_rasters.csv` allows additional Earth Engine assets; the committed fruit-tree-suitability row is disabled.

**Prediction grid.** Cells are sampled at approximately 10 km (`PREDICTION_GRID_SCALE_M = 10000`) on land inside the study area. Prediction years in Python are 2020–2025. Intended merged tables are `data/dataset2.csv`, `data/prediction_grid_10km.csv`, and `data/prediction_grid_covariates_2020_2025.csv` (all gitignored).

### 2.6 Boosted regression trees (script 03)

Script `03_train_predict_brt_simple.R` reads `dataset2.csv` and the prediction-grid covariate table. Hansen land missing values are filled with 0. It builds `N_DATASETS = 50` resampled tables, each containing all events and `CONTROLS_PER_EVENT = 50` controls sampled without replacement (`RANDOM_SEED = 20260826`). Each table is passed to `dismo::gbm.step` with Bernoulli family, tree complexity 5, learning rate 0.0012, bag fraction 0.7, and 10 internal folds [15]. Ensemble predictions on the grid are summarised as min, max, mean, and median probability, then rasterised (`terra`, EPSG:4326) to annual GeoTIFFs and PNG maps.

Apparent (resubstitution) performance on 2020–2025 training points is optional (`RUN_PERFORMANCE_EVALUATION <- TRUE`) and uses `pROC`, `PRROC`, and `caret`. The default threshold maximises the product of sensitivity and positive predictive value. The script labels these metrics “apparent”; they are not temporal hold-outs. Optional TreeSHAP (`treeshap`) is computed on the GBM link scale and averaged across selected ensemble members [23]. Relative-influence tables and partial-dependence (marginal-effect) curves are written under `outputs/model_diagnostics/`.

Software named for this script includes R packages `gbm`, `dismo`, `terra`, `sf`, `ggplot2`, `pROC`, `PRROC`, `caret`, and `treeshap`.

### 2.7 Temporal forward validation (script 03b)

Script `03b_temporal_forward_validation_brt.R` is designed to emulate operational use: train only on years before each target year, then predict the held-forward year. Documented folds are train 2001–2021 → predict 2022; 2001–2022 → 2023; 2001–2023 → 2024; and 2001–2024 → 2025 (`TEMPORAL_TARGET_YEARS <- 2022:2025`). Hyperparameters match script 03 except that each fold uses 10 resampled datasets and 100 controls per event. Outputs are written under `models/tfv/` and `outputs/tfv/` so they do not overwrite the all-years run. Performance CSVs and ROC/PR figures are specified but not deposited.

### 2.8 SuperLearner prototype (script 05)

Script `05_train_predict_SuperLearner_CV.R` is labelled a first-pass prototype. It fits a SuperLearner [24] with library `SL.glm`, `SL.rpart_tuned`, and `SL.ranger_tuned`, internal 5-fold stacking (`method.NNloglik`), and 10 outer folds. Three tuning candidates (`quick_conservative`, `quick_balanced`, `quick_flexible`) vary rpart complexity and ranger `mtry`/trees. A threshold grid from 0.001 to 0.50 is scored by F1. Hansen missing values are zero-filled; remaining missing predictors are median-imputed. Prediction years are again 2020–2025. Additional packages: `SuperLearner`, `rpart`, `ranger`.

### 2.9 Interactive reports

`04_interactive_prediction_report_Shiny.Rmd` is a Shiny/leaflet viewer for annual rasters. `04b_interactive_prediction_report_email.Rmd` knits a static HTML report. Both expect local `outputs/` rasters that Git does not contain.

### 2.10 What this draft computed

We cloned `https://github.com/BoazBaliejukia/KSPH_Code` on 3 September 2026 and read every tracked file under `R_python_code/` and `config/`. Descriptive counts in the Results are from `config/outcomes_csv.csv` only. We did not run `01_make_dataset1.R` (that would create stochastic pseudo-absences, not new scientific measurements of spillover), did not call Earth Engine, and did not fit BRT or SuperLearner models. Ethics: the draft uses only publicly posted village- or health-zone-level coordinates and script metadata. No new human-subjects data were collected.

---

## 3. Results

### 3.1 Composition of the committed presence linelist

`config/outcomes_csv.csv` has 27 data rows and six columns (`ID`, `Type`, `Country`, `Year`, `Lat`, `Long`). Type codes are Z (n = 18), S (n = 6), and B (n = 3). Country codes are DRC (11), Republic of the Congo (6), Uganda (6), Gabon (2), South Sudan (1), and Guinea (1). Years with at least one record are 2001–2005, 2007–2008, 2011–2014, 2017–2018, 2020, 2022, 2025, and 2026. Three events fall in 2012 (the modal year); two events fall in each of 2001, 2002, 2003, 2004, 2007, 2018, 2022, and 2025.

Coordinates span 5.25956°S to 8.6225°N and 10.06423°W to 32.72762°E. The westernmost point is Meliandou, Guinea (2013, Z). The easternmost is Nakisamata, Uganda (2011, S). The southernmost points are Bamoukamba (2007) and Kaluamba (2008) in the DRC. All 26 rows with year 2001–2025 lie inside `STUDY_AREA_BBOX`.

**Table 1.** Presence records in `config/outcomes_csv.csv`.

| ID | Type | Country | Year | Latitude | Longitude |
|---|---|---|---:|---:|---:|
| GEB2 | Z | Gabon | 2001 | 0.05503 | 11.81095 |
| Medemba | Z | Gabon | 2001 | 0.70011 | 14.15803 |
| Entsiami | Z | ROC | 2002 | 0.13332 | 14.26667 |
| Lossi1 | Z | ROC | 2002 | 0.29870 | 14.50750 |
| Mvoula | Z | ROC | 2003 | 0.06823 | 14.41997 |
| Mbanza2 | Z | ROC | 2003 | 0.56186 | 14.65726 |
| Yambio | S | SSudan | 2004 | 4.55676 | 28.40155 |
| Lokoue | Z | ROC | 2004 | 0.90638 | 15.17512 |
| Etoumbi | Z | ROC | 2005 | 0.49444 | 14.67861 |
| Bamoukamba | Z | DRC | 2007 | −5.25956 | 21.40954 |
| Bundibugyo | B | Uganda | 2007 | 0.70383 | 30.11752 |
| Kaluamba | Z | DRC | 2008 | −5.24009 | 21.41030 |
| Nakisamata | S | Uganda | 2011 | 0.64443 | 32.72762 |
| Kibaale | S | Uganda | 2012 | 0.62140 | 31.16849 |
| Isiro | B | DRC | 2012 | 2.77176 | 27.61961 |
| Luwero | S | Uganda | 2012 | 0.57843 | 32.54796 |
| Meliandou | Z | Guinea | 2013 | 8.62250 | −10.06423 |
| Inkanamongo | Z | DRC | 2014 | −0.71387 | 20.53024 |
| Likati | Z | DRC | 2017 | 3.29898 | 23.54304 |
| Bikoro | Z | DRC | 2018 | −0.73684 | 18.42141 |
| NorthKivu | Z | DRC | 2018 | 0.60589 | 29.30645 |
| Mbandaka | Z | DRC | 2020 | 0.02998 | 18.28004 |
| Wangata | Z | DRC | 2022 | 0.04135 | 18.27702 |
| Madudu | S | Uganda | 2022 | 0.66965 | 31.46862 |
| Bulape | Z | DRC | 2025 | −4.67627 | 21.60698 |
| Wakiso | S | Uganda | 2025 | 0.49298 | 32.50361 |
| Bunia | B | DRC | 2026 | 1.61247 | 30.25402 |

Rows through 2022 match the village identifiers and coordinates used in the Telford et al. replication table [1,14]. The three later rows are new relative to that 2001–2022 archive.

### 3.2 Year filter and the 2026 record

`01_make_dataset1.R` keeps presence years only if they belong to `ABSENCE_YEARS` (`2001:2025`). The Bunia 2026 row is therefore **dropped before pseudo-absence binding**. After that filter the training-eligible linelist has **26 events** (18 Z, 6 S, 2 B). The 2026 Bundibugyo outbreak is visible in the raw CSV and in WHO situation reports [20,21], but it is not in the script’s current modelling window. Temporal forward validation similarly ends at target year 2025. Anyone wishing to score the Ituri outbreak on a predicted surface would need to extend `ABSENCE_YEARS`, the prediction-year vectors, and the Earth Engine export years; those changes are not in the committed defaults.

### 3.3 What the three post-2022 points represent, relative to outbreak reports

We did not treat the new coordinates as independently verified index-case villages.

- **Bulape, 2025 (Z, DRC).** WHO reported that the 2025 Kasai outbreak was confined to six health areas in Bulape Health Zone [18,19]. The linelist coordinate (−4.676, 21.607) is consistent with that health zone. Previous Kasai outbreaks in 2007 and 2008 appear in the same file as Bamoukamba and Kaluamba, a few tens of kilometres away [18].

- **Wakiso, 2025 (S, Uganda).** WHO listed confirmed or probable Sudan virus disease in Fort Portal City, Jinja, Kampala, Kyegegwa, Mbale, Ntoroko, and Wakiso [17]. Public accounts of the index patient describe care-seeking across Mbale, Kampala, and Wakiso [17]. The single Wakiso coordinate is therefore a district-level placeholder, not a demonstration that spillover occurred at that point.

- **Bunia, 2026 (B, DRC).** WHO’s initial Disease Outbreak News stated that the event was thought to have originated in Mongbwalu Health Zone, with subsequent cases in Rwampara and Bunia [20]. Later updates described Bunia as one of the most affected health zones in a multi-province epidemic [21]. The linelist stores one 2026 point in Bunia. That is a reasonable label for “Ituri / Bunia-area outbreak” and is **not** a published index-household coordinate.

These distinctions matter because the BRT outcome is the location at which an event is *recorded*, which may be a village of residence, a health zone of care, or an outbreak centroid [1].

### 3.4 Intended but undeposited analysis products

The public tree tracks `R_python_code/`, `config/`, `README.md`, and `.gitignore`. The ignore file excludes `data/`, `models/`, and `outputs/` by design. Consequently, the following objects **do not exist in the GitHub snapshot** and were not generated here:

1. `data/dataset1.csv` — designed to contain 26 presence rows plus 5,000 pseudo-absences after the year filter.
2. `data/dataset2.csv` and `data/prediction_grid_covariates_2020_2025.csv` — designed to hold the 35 model predictors at training points and on the 10 km grid.
3. `models/model_list.rds` and `models/dsl.rds` — the 50 BRT fits and their resampled training tables.
4. Annual probability, ROR, one-year-change, and top-1% rasters for 2020–2025.
5. Apparent ROC/PR tables and caret confusion matrices under `outputs/performance_evaluation/`.
6. Temporal-forward-validation metrics for 2022–2025 under `outputs/tfv/`.
7. SuperLearner cross-validated predictions and 2020–2025 maps under `outputs/predictions/superlearner_cv/`.
8. TreeSHAP tables and Shiny/HTML reports.

We therefore report **no** AUC, sensitivity, specificity, F1, relative-influence ranking, partial-dependence shape, or mapped ROR value. Those quantities are what the scripts are written to produce after a local Earth Engine export. Stating any of them as results of this run would be fabrication.

### 3.5 Comparison of pipeline defaults with Telford et al. (2025)

Table 2 records committed differences that can be read from scripts without rerunning models.

**Table 2.** Selected design choices in KSPH Code versus the published 2025 analysis [1,14].

| Choice | Telford et al. (2025) / AnnualEbolaPrediction_EID | KSPH Code (this repository) |
|---|---|---|
| Presence years in public table | 24 events, 2001–2022 | 27 events, 2001–2026 (26 enter `dataset1`) |
| Background sampling | 10,000 LandScan-weighted points | 5,000 uniform land points |
| Forest product | Hansen GFC 2023 v1.11 in the replication scripts | Hansen GFC 2025 v1.13 |
| Climate | FLDAS precipitation seasonality; MODIS LST | ERA5-Land precipitation, temperature, PET; MODIS NDVI |
| Spatial rings | 10, 25, 50, 100, 150 km | 10, 25, 50, and labelled 100 km (≈75 km outer) |
| Ensemble | 100 BRTs; 50 absences per presence; complexity 3; learning rate 0.0009–0.001 | 50 BRTs; 50 controls per event; complexity 5; learning rate 0.0012 |
| Validation in public scripts | Leave-one-year-out over 2001–2021 | Apparent 2020–2025 metrics plus temporal forward folds 2022–2025 |
| Additional learner | Not in the public EID scripts | SuperLearner prototype (glm / rpart / ranger) |
| Prediction years in scripts | 2021–2022 grids | 2020–2025 grids |
| Deposited covariates | Yes (`tdf_clean.csv`) | No (exports gitignored) |

The 2022 Wangata (DRC, Z) and Madudu/Mubende (Uganda, S) rows that Telford et al. held out of training [1] are present in the new linelist and would enter an all-years KSPH fit unless a user sets a type or year filter.

---

## 4. Discussion

### 4.1 A pipeline that can be updated after 2022, and why that matters

Telford et al. showed that annually refreshed forest-loss, fragmentation, and population layers can rank equatorial African cells by relative spillover odds, and that the two 2022 events sat in the extreme tail of those ranks or of 2021–2022 change [1]. The ecological rationale—forest edge, recent loss, and human presence as correlates of where index cases are found—is consistent with earlier niche and fragmentation studies [8–11]. The scientific gap after 2025 is not another static suitability map. It is whether a documented, re-runnable pipeline can ingest new events and new annual rasters without rewriting the analysis from scratch.

KSPH Code is a direct attempt at that gap. It separates linelist construction, Earth Engine extraction, model fitting, temporal hold-forward testing, and reporting. It upgrades Hansen to the 2025 v1.13 release, replaces FLDAS/LST seasonality with ERA5-Land totals and anomalies plus NDVI [22], and adds an explicit 2022–2025 forward-validation loop. Those are methodological improvements on paper. They are not yet empirical improvements: without `dataset2.csv` and `model_list.rds`, we cannot say whether ERA5 anomalies or the SuperLearner library outperform the published BRT ensemble, nor whether Bulape 2025 or Wakiso 2025 would have ranked in the top 0.1% of a 2024-trained ROR surface—the analogue of the 2022 assessment in Telford et al. [1].

### 4.2 The 2025–2026 outbreaks as a test of honesty, not as model output

The new linelist rows are epidemiologically real events [17–21]. Using them as presence points is reasonable for an updated niche model, but they are not interchangeable with the 2001–2021 village series. The 2025 Uganda SVD outbreak was a multi-district chain after a mobile index patient [17]. The 2025 Kasai outbreak was spatially compact in Bulape Health Zone and was contained within three months [18,19]. The 2026 Bundibugyo outbreak, by contrast, became a large, multi-province epidemic without a licensed vaccine [20,21]. A presence–background model that treats each outbreak as a single point-year cannot represent spatial expansion or secondary transmission. It can, at best, ask whether the *recorded origin or label location* sat in a landscape that the previous year’s covariates would have ranked highly.

That question is precisely what script 03b is written to answer for 2022–2025. Until those folds are run, the public-health implication of this repository is limited: the linelist is ready; the scoring is not. Claiming that the pipeline “predicted” Bulape or Wakiso would be false.

### 4.3 Presence–background design and causal reading

Random land pseudo-absences differ from the LandScan-weighted background in Telford et al. [1,13]. Uniform-space controls make forested, sparsely populated cells look more distinctive relative to background than they would against a population-weighted null; population effects will be estimated on a different contrast. Neither design identifies a causal effect of deforestation. Mutually adjusted landscape terms in a BRT or SuperLearner should not be read as independent etiologic effects (the “Table 2 fallacy”) [16]. Predictions refer to where an event is recorded, which may be residence, care-seeking, or a health-zone centroid [1]. Undetected spillovers remain coded as absences.

BRT ensembles remain an appropriate class for this prediction problem because they handle nonlinearities and correlated environmental layers [15]. SuperLearner is a theoretically grounded way to stack GLM, trees, and random forests [24], but the committed library is a “quick” three-learner grid and is described in the script header as a first pass. Apparent AUC on the same years used to train the ensemble (script 03) will overstate operational accuracy; the temporal script exists because of that problem.

### 4.4 Limitations of this draft and of the public archive

The presence sample is 27 rows (26 after the year filter). That is enough to document a linelist and a recipe; it is not enough to support fine-grained species-specific maps without extreme uncertainty, even after models are fit. Judson et al. showed that ecological context differs among Ebolavirus species [4]; the default KSPH run pools Z, S, and B. The 2026 Bunia point is epidemiologically important and currently unused. LandScan is frozen at 2023, so 2024–2025 population layers are last-available-year copies. The `50_100km` ring label does not match the 0.674° outer radius. Prediction-grid extraction omits the outermost training ring, so training and grid predictor sets can diverge unless the notebook is aligned. Forest-loss lags are structurally missing in 2001–2002. We did not audit buffer geometry on Earth Engine.

Most importantly, this draft cannot report discrimination, calibration, variable importance, or maps. Those are the Results a journal would expect from a completed niche-model paper. Their absence is a property of the public snapshot, not an omission we filled with invented numbers.

---

## 5. Conclusion

KSPH Code specifies a complete, annually updateable event-niche pipeline for equatorial Africa: a 2001–2026 Ebolavirus presence linelist, land-only pseudo-absence sampling, Hansen/LandScan/ERA5-Land/MODIS/SRTM extraction over multi-scale donut buffers, a 50-member Bernoulli BRT ensemble with ROR mapping for 2020–2025, a 2022–2025 temporal-forward-validation loop, and a SuperLearner prototype. The committed linelist adds three outbreaks after Telford et al. (2025)—Bulape (Zaire, DRC, 2025), Wakiso (Sudan, Uganda, 2025), and Bunia (Bundibugyo, DRC, 2026)—of which the first two enter the current 2001–2025 training window. Covariate tables, fitted models, validation metrics, and rasters are intentionally local and are not in Git. Until those artefacts are produced or deposited, the honest scientific product of the public repository is a methods specification and an updated presence file, not a new estimate of where spillover is most likely.

---

## Data and code availability

- Upstream repository: https://github.com/carsontelford/KSPH-Code  
- Fork used for this draft (3 September 2026 snapshot): https://github.com/BoazBaliejukia/KSPH_Code  
- Presence linelist: `config/outcomes_csv.csv` (27 rows).  
- Predictor list: `config/predictor_list.csv` (35 model covariates).  
- Study-area geometry: `config/africacountries_nolakes.shp`.  
- Scripts: `R_python_code/01_make_dataset1.R`, `02_predGrid_trainSet_extraction.ipynb`, `02_drivehelper_covariates.py`, `02_prediction_grid_helper.py`, `03_train_predict_brt_simple.R`, `03b_temporal_forward_validation_brt.R`, `04_interactive_prediction_report_Shiny.Rmd`, `04b_interactive_prediction_report_email.Rmd`, `05_train_predict_SuperLearner_CV.R`.  
- Related published analysis: Telford et al., *Emerg Infect Dis.* 2025;31(4):689–698. doi:[10.3201/eid3104.241193](https://doi.org/10.3201/eid3104.241193) [1].

Generated `data/`, `models/`, and `outputs/` directories are **not** in the public repository.

---

## References

1. Telford CT, Amman BR, Towner JS, Montgomery JM, Lessler J, Shoemaker T. Predictive model for estimating annual Ebolavirus spillover potential. *Emerg Infect Dis.* 2025;31(4):689-698. doi:10.3201/eid3104.241193

2. Rouquet P, Froment JM, Bermejo M, Kilbourn A, Karesh W, Reed P, et al. Wild animal mortality monitoring and human Ebola outbreaks, Gabon and Republic of Congo, 2001–2003. *Emerg Infect Dis.* 2005;11(2):283-290. doi:10.3201/eid1102.040533

3. Leroy EM, Kumulungui B, Pourrut X, Rouquet P, Hassanin A, Yaba P, et al. Fruit bats as reservoirs of Ebola virus. *Nature.* 2005;438(7068):575-576. doi:10.1038/438575a

4. Judson SD, Fischer R, Judson A, Munster VJ. Ecological contexts of index cases and spillover events of different ebolaviruses. *PLoS Pathog.* 2016;12(8):e1005780. doi:10.1371/journal.ppat.1005780

5. Leroy EM, Epelboin A, Mondonge V, Pourrut X, Gonzalez JP, Muyembe-Tamfum JJ, et al. Human Ebola outbreak resulting from direct exposure to fruit bats in Luebo, Democratic Republic of Congo, 2007. *Vector Borne Zoonotic Dis.* 2009;9(6):723-728. doi:10.1089/vbz.2008.0167

6. Centers for Disease Control and Prevention. History of Ebola disease outbreaks. Accessed as cited by Telford et al. [1]. https://www.cdc.gov/vhf/ebola/history/chronology.html

7. Mylne A, Brady OJ, Huang Z, Pigott DM, Golding N, Kraemer MUG, et al. A comprehensive database of the geographic spread of past human Ebola outbreaks. *Sci Data.* 2014;1:140042. doi:10.1038/sdata.2014.42

8. Pigott DM, Millear AI, Earl L, Morozoff C, Han BA, Shearer FM, et al. Updates to the zoonotic niche map of Ebola virus disease in Africa. *eLife.* 2016;5:e16412. doi:10.7554/eLife.16412

9. Schmidt JP, Park AW, Kramer AM, Han BA, Alexander LW, Drake JM. Spatiotemporal fluctuations and triggers of Ebola virus spillover. *Emerg Infect Dis.* 2017;23(3):415-422. doi:10.3201/eid2303.160101

10. Rulli MC, Santini M, Hayman DTS, D’Odorico P. The nexus between forest fragmentation in Africa and Ebola virus disease outbreaks. *Sci Rep.* 2017;7:41613. doi:10.1038/srep41613

11. Olivero J, Fa JE, Real R, Márquez AL, Farfán MA, Vargas JM, et al. Recent loss of closed forests is associated with Ebola virus disease outbreaks. *Sci Rep.* 2017;7:14291. doi:10.1038/s41598-017-14727-9

12. Hansen MC, Potapov PV, Moore R, Hancher M, Turubanova SA, Tyukavina A, et al. High-resolution global maps of 21st-century forest cover change. *Science.* 2013;342(6160):850-853. doi:10.1126/science.1244693

13. Sims K, Reith A, Bright E, Kaufman J, Pyle J, Epting J, et al. LandScan Global 2022. Oak Ridge (TN): Oak Ridge National Laboratory; 2023.

14. Telford CT. AnnualEbolaPrediction_EID. GitHub. 2024. https://github.com/carsontelford/AnnualEbolaPrediction_EID

15. Elith J, Leathwick JR, Hastie T. A working guide to boosted regression trees. *J Anim Ecol.* 2008;77(4):802-813. doi:10.1111/j.1365-2656.2008.01390.x

16. Westreich D, Greenland S. The table 2 fallacy: presenting and interpreting confounder and modifier coefficients. *Am J Epidemiol.* 2013;177(4):292-298. doi:10.1093/aje/kws412

17. World Health Organization. Sudan virus disease – Uganda. Disease Outbreak News. 26 April 2025. https://www.who.int/emergencies/disease-outbreak-news/item/2025-DON566

18. World Health Organization. Ebola virus disease – Democratic Republic of the Congo. Disease Outbreak News. 1 December 2025. https://www.who.int/emergencies/disease-outbreak-news/item/2025-DON589

19. Journal of Epidemiology and Global Health. Ebola in Kasai, Democratic Republic of the Congo: lessons from a contained outbreak for epidemic preparedness and global health security. *J Epidemiol Glob Health.* 2026. doi:10.1007/s44197-026-00616-w

20. World Health Organization. Ebola disease caused by Bundibugyo virus, Democratic Republic of the Congo & Uganda. Disease Outbreak News. 2026. https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON602

21. World Health Organization. Ebola disease caused by Bundibugyo virus – Democratic Republic of the Congo. Disease Outbreak News. 2026. https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON614

22. Muñoz-Sabater J, Dutra E, Agustí-Panareda A, Albergel C, Arduini G, Balsamo G, et al. ERA5-Land: a state-of-the-art global reanalysis dataset for land applications. *Earth Syst Sci Data.* 2021;13:4349-4383. doi:10.5194/essd-13-4349-2021

23. Telford CT. KSPH-Code: event ecological niche pipeline. GitHub. 2026. https://github.com/carsontelford/KSPH-Code

24. van der Laan MJ, Polley EC, Hubbard AE. Super Learner. *Stat Appl Genet Mol Biol.* 2007;6:Article 25. doi:10.2202/1544-6115.1309
