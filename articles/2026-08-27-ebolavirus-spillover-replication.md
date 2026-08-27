# Landscape covariates at Ebolavirus spillover sites in equatorial Africa, 2001–2022: pipeline documentation and a descriptive re-analysis of open replication data

**Manuscript type.** Methods and data note (unpublished draft).  
**Provenance.** This draft documents the public replication repository for Telford et al. (2025), *Emerging Infectious Diseases* [1], as archived at [carsontelford/AnnualEbolaPrediction_EID](https://github.com/carsontelford/AnnualEbolaPrediction_EID) and mirrored at [BoazBaliejukia/AnnualEbolaPrediction_EID](https://github.com/BoazBaliejukia/AnnualEbolaPrediction_EID). The original boosted-regression-tree (BRT) ensemble, leave-one-year-out cross-validation metrics, and 2021–2022 relative-odds maps were produced by Telford and colleagues. The descriptive tables in the Results are independently computed from the deposited file `tdf_clean.csv`. Fitted ensemble objects and prediction-grid rasters are not in the public repository and were **not** regenerated for this draft.

---

## Abstract

**Background.** Zoonotic spillover of Ebolavirus species into humans remains infrequent but consequential. Forest cover, forest loss and fragmentation, human population distribution, and meteorologic seasonality have been linked to where index cases are identified. Telford et al. (2025) fitted annually updated BRT ensembles to estimate relative odds of spillover across equatorial Africa.

**Objective.** To document the open analysis pipeline in sufficient detail for independent use, and to describe the deposited presence–absence table without fabricating model-output maps or performance statistics that cannot be recomputed from public files.

**Methods.** We reviewed the three R scripts (absence-point generation; model fitting and cross-validation; prediction over a 1 km grid), Google Earth Engine extraction scripts, and the post-processed table `tdf_clean.csv` (24 documented spillover events during 2001–2022 and 10,000 population-weighted absence locations). Following the published training window, we restricted independent descriptives to 2001–2021 (22 presence rows). Covariates were reconstructed with the same derived forest-loss, log-population, fragmentation, and scaling steps as in the fitting script. Presence and absence locations were compared with medians, interquartile ranges, and unadjusted two-sided Mann–Whitney tests.

**Results.** Same-year forest cover within 10 km was higher at spillover villages (median 79.2%, IQR 53.7–97.2) than at absence locations (median 35.7%, IQR 16.4–59.6). Same-year forest-loss and fragmentation scores were likewise higher at presence sites, and precipitation seasonality was lower. Log human population within 10 km did not differ detectably, whereas log population within 150 km was lower at presence sites. Two 2022 events (coordinates consistent with the Wangata, Democratic Republic of the Congo, and Mubende, Uganda, outbreaks) are present in the table and were excluded from the training filter `Year > 2000 & Year < 2022`. The public scripts document a leave-one-year-out ensemble AUC of 0.88 and a 2×2 classification of 19 true positives, 3 false negatives, 2,731 false positives, and 7,269 true negatives among 22 training events and 10,000 absences; those figures match the published multispecies sensitivity (86.4%) and specificity (72.7%). We did not refit the 100-member BRT ensemble or the 2021–2022 prediction grids, because the required raster stacks and `.RData` model objects are not deposited.

**Conclusions.** The open table supports a transparent description of landscape context at documented spillover villages relative to population-weighted background points. It does not, by itself, reproduce annual relative-odds maps. Independent replication of the predictive product still requires the large covariate rasters used in scripts P1 and P3.

**Keywords.** Ebolavirus; zoonotic spillover; forest fragmentation; boosted regression trees; equatorial Africa; open replication.

---

## 1. Introduction

Outbreaks of Ebola virus disease (EVD) typically begin when an Ebolavirus species infects a person from an animal source and then spreads in human networks [1,2]. Four species have caused human disease: *Zaire ebolavirus* (EBOV, often abbreviated ZEBOV in the source paper), *Sudan ebolavirus* (SUDV), *Bundibugyo ebolavirus* (BDBV), and *Taï Forest ebolavirus* (TAFV) [1]. A definitive reservoir has not been established for each species, although forest-associated bats are leading candidates and experimental and outbreak investigations have implicated fruit bats in at least some transmission settings [3–5]. Index-case investigations have also implicated hunting, butchering, and market handling of wildlife [2,5].

The public-health burden of secondary transmission can be extreme. The 2014–2016 West African epidemic caused nearly 30,000 reported cases and more than 11,000 deaths, and the 2018–2020 outbreak in the Democratic Republic of the Congo (DRC) caused approximately 3,500 cases and 2,300 deaths [1]. Identification of spillover events has become more frequent in the twenty-first century than in the first two decades after 1976, reflecting some combination of improved detection and changing human–animal contact [1,6]. Because spillover remains rare in any given place and year, surveillance cannot be equally intensive everywhere. Models that rank locations by relative odds of an index case, using covariates that are updated annually, are therefore of operational interest [1,7].

Several ecological studies have associated EVD emergence with forest structure and climate. Pigott and colleagues mapped a zoonotic niche for Ebola virus disease using environmental suitability [8]. Schmidt et al. related human and great-ape spillover intensity to vegetation, population density, and rainfall seasonality, reporting higher intensity at wet–dry transitions and at both very low and very high human population densities [9]. Rulli et al. found that index cases during 2004–2014 occurred preferentially in fragmentation hotspots rather than in hotspots of forest *loss* per se [10]. Olivero et al. reported that outbreaks along rainforest margins were associated with closed-forest loss in the preceding two years [11]. Those results motivate models that treat forest loss and fragmentation as time-varying predictors, not only as static habitat layers.

Telford et al. addressed a practical gap: high-resolution forest-change and population datasets are released annually, typically with a lag of several months, and can be used to refresh relative-odds surfaces rather than to produce a single historical suitability map [1,12,13]. They fitted BRT ensembles for 2001–2021 and generated relative odds ratios (RORs) for 2021 and 2022, holding out the two 2022 spillover events for assessment [1]. Code and a cleaned presence–absence table were deposited on GitHub [14].

BRT methods combine regression trees with boosting and are widely used for ecological prediction when many correlated predictors and nonlinear effects are expected [15]. They are not a substitute for etiologic identification of a reservoir or a causal decomposition of land-use effects [1,16]. The present draft has two aims: (i) to write Methods from the actual scripts and data files, and (ii) to report independently computed descriptives from `tdf_clean.csv`, while clearly separating those descriptives from published ensemble metrics and maps that we did not regenerate.

---

## 2. Methods

### 2.1 Study design and outcome

The analysis implemented in the repository is a presence–background (presence–absence) predictive study of documented Ebolavirus spillover into humans. The outcome in script `P2 Model fit and CV.R` is a binary indicator `spillover`, coded 1 when `Strain` is Zaire, Sudan, or Bundibugyo and 0 when `Strain` is `Absence_logpdw`. TAFV is excluded because the only documented spillover was in 1994, before the forest-change and population series used here [1].

Presence points are geographic coordinates of the village of residence of outbreak index cases that were not attributed to recrudescence from persistently infected survivors [1]. The README states that the deposited table contains 24 spillover events from 2001–2022 and 10,000 random absence locations [14]. Script P2 restricts model training to `Year > 2000 & Year < 2022` (calendar years 2001–2021), leaving 2022 events out of the fitting window.

The intended predictive target in script `P3 Prediction grid.R` is not a probability calibrated to absolute annual incidence, but a relative odds surface: fitted odds in each 1 km cell divided by the mean fitted odds across the study area in that year (ROR), as defined in the published paper [1].

### 2.2 Study area

Predictors were compiled for continental Africa within 10 degrees of the Equator, excluding locations with historical mean annual precipitation below 500 mm, on the rationale that documented spillover has occurred in vegetation-dense settings [1,9]. Absence-point generation in `P1 Absence point generation.R` samples from a LandScan-derived population raster covering that domain. In `tdf_clean.csv`, background coordinates span approximately 10°S–10°N and 14°W–40°E.

### 2.3 Data sources

**Spillover linelist.** Script P1 reads `EVD_spill_csv.csv` (not deposited) and keeps records with `Include != "No"`, using confirmed longitude and latitude. Outbreaks without published coordinates were resolved with CDC outbreak-response consultation in the original study [1].

**Forest cover and loss.** Google Earth Engine scripts load Hansen Global Forest Change (`UMD/hansen/global_forest_change_2023_v1_11`) and compute mean tree cover within donut buffers of 10, 25, 50, 100, and 150 km [12]. Same-year cover (`fcsy*`), one-year-prior cover (`fc1yp*`), two-year-prior cover (`fc2yp*`), and three-year-prior cover (`fc3yp*`) are stored in `tdf_clean.csv`. Script P2 defines forest-loss layers as differences between successive cover layers (for example `flsy10km = fc1yp10km - fcsy10km`), then multiplies loss by 10,000 and cover by 100 for numerical scaling.

**Fragmentation.** Forests were classified using a ≥70% cover threshold (sensitivity analysis at 80% in the published paper) and then labelled as fragmented following the patch, transitional, edge, and perforated typology of Riitters et al. [1,17]. Columns `frag.7allarea*` in the CSV are multiplied by 100 in P2.

**Human population.** Engine scripts sum Oak Ridge National Laboratory LandScan Global counts within the same buffers [13]. P2 uses `log(pcsy* + 1)` at each radius and a product term `logpcsy10km * fcsy10km` (`logpcfc10km`) to represent co-occurrence of people and remaining forest.

**Climate and terrain.** Potential evapotranspiration (`pet`), night-time land-surface temperature (`ntlst`), elevation, temperature seasonality, and precipitation seasonality are stored as 10 km means. A MODIS MOD11A2 script extracts land-surface temperature; an FLDAS script extracts precipitation seasonality. P2 divides `pet` by 1,000 and `tempseason` by 10.

P1 and P3 state that the raw rasters and prediction-grid extracts are too large to host on GitHub. Only `tdf_clean.csv` is provided as the analysis-ready table for script P2 [14].

### 2.4 Absence locations

P1 generates 10,000 background points with `seegSDM::bgSample` on `log(population_density + 1) + 1`, with `prob = TRUE` and `replace = FALSE`, so sampling is weighted toward populated cells. Each point is assigned a year drawn uniformly from 2001–2021 and labelled `Absence_logpdw`. The published paper describes this weighting as an adjustment for higher reporting probability in populated areas; a sensitivity analysis used health-facility coordinates as alternative background points [1,18]. Facility-based absences are not in `tdf_clean.csv`.

### 2.5 Predictive models in the source scripts

After deriving scaled covariates, P2 retains 36 predictors (`gbm.x = 7:42`) plus the binary outcome (`gbm.y = 43`) and coordinates (`gbm.coords = 3:4`). The published paper states that 30 of these covariates are spatially scaled forest and population layers and six are unscaled 10 km layers [1].

**Leave-one-year-out cross-validation (documented, not rerun here).** For each calendar year used as a test fold, training data exclude that year. An ensemble of 100 datasets is built by sampling 50 absence rows per remaining presence row (`presnum * 50`, with replacement in the CV loop). Each dataset is passed to `gbm.step` (Bernoulli BRT; tree complexity 3; learning rate 0.001; bag fraction 0.5) in parallel via `snowfall` (`sfLapply`, 4 CPUs) and `seegSDM` [15,19]. Test-fold predictions are averaged across the 100 models. Concatenated out-of-year fitted values are evaluated with `pROC::roc`. Sensitivity and specificity use a cutoff that maximises the product of the two, via `ROCR` [1].

**Final ensemble (documented, not rerun here).** The full 2001–2021 table is resampled similarly (`replace = FALSE` in the final-fit loop) and fitted with `seegSDM::runBRT` (learning rate 0.0009; 3-fold internal CV). The script writes `FittedModelEnsembles/AllspeciesFullModel_WeightedAbsence.RData` and later `load()`s a local CDC path that is not in the repository. Relative-influence medians across the 100 models and partial-dependence curves are computed from those objects.

**Prediction (documented, not rerun here).** P3 joins Google Earth Engine extracts for 2021 and 2022 grids, rebuilds the same derived covariates, predicts from each ensemble member, and aggregates to cell-level odds and RORs. Those grid CSVs are not deposited.

Software named in the scripts includes R with `raster`, `tidyverse`, `sf`, `sp`, `rgdal`, `tmap`, `dismo`, `exactextractr`, `rgeos`, `gbm`, `seegSDM`, `maptools`, `fields`, `parallel`, `snowfall`, `pROC`, `ROCR`, `gridExtra`, and `grid`. Several of those packages (`rgdal`, `rgeos`, `maptools`) have been retired from CRAN, which is a practical barrier to bit-for-bit reruns.

### 2.6 Independent descriptive analysis for this draft

We read `tdf_clean.csv` from the 27 August 2026 snapshot of the BoazBaliejukia fork (identical to the upstream Telford repository). Derived variables were created with the same formulae and scalings as P2. Training-window rows were `Year` 2001–2021. For each covariate we report presence and absence sample size, median, and 25th–75th percentiles. Two-sided Mann–Whitney *U* tests are exploratory and unadjusted for the number of covariates; they are not causal estimates and are not a reproduction of BRT relative influence [16]. Analyses used Python 3.12 with pandas 2.x and SciPy.

We did not fit BRT, logistic, or other predictive models for this draft, because a simplified refit would not be the published 100-member ensemble and could be mistaken for a replication of Telford et al.

### 2.7 Ethics

This draft uses only publicly posted, aggregated geographic coordinates and remotely sensed covariates. No new human-subjects data were collected.

---

## 3. Results

### 3.1 Composition of `tdf_clean.csv`

The deposited table has 10,024 rows and 41 columns. Strain codes are `Absence_logpdw` (n = 10,000), `Zaire` (n = 17), `Sudan` (n = 5), and `Bundibugyo` (n = 2), matching the README total of 24 presence events [14]. Absence years run from 2001 through 2021 only (445–523 absences per year). Two presence rows are dated 2022 and are dropped by the P2 training filter, leaving **22 presence and 10,000 absence rows** for 2001–2021 (16 Zaire, 4 Sudan, 2 Bundibugyo).

Presence coordinates in the training window range from 5.26°S to 8.62°N and from 10.06°W to 32.73°E. The westernmost training event is a 2013 Zaire-coded point at 10.06°W, 8.62°N (West Africa). Eastern training events include Sudan-coded points in 2004, 2011, and 2012 beyond 28°E.

The two 2022 hold-out rows are:

| Year | Strain | Longitude (`x`) | Latitude (`y`) |
|------|--------|-----------------|----------------|
| 2022 | Zaire  | 18.27702 | 0.04135 |
| 2022 | Sudan  | 31.46862 | 0.66965 |

These coordinates lie in the regions of the 2022 Wangata (Mbandaka, DRC) EBOV outbreak and the Mubende (Uganda) SUDV outbreak described by Telford et al. [1]. We did not independently re-verify village names against outbreak reports beyond the coordinate match to that description.

Lagged forest-cover fields are incomplete early in the series: `fc2yp*` is missing for 496 rows (all year 2001) and `fc3yp*` for 1,004 rows (years 2001–2002). Consequently, complete cases on the full 36-covariate set number 9,018 of 10,022 training-window rows (18 of 22 presence rows; 9,000 of 10,000 absences). Same-year cover, population, fragmentation, and climate fields have no missing values in the CSV.

### 3.2 Landscape covariates at presence versus absence locations

Table 1 summarises selected covariates after P2 scaling. Same-year forest cover was substantially higher at spillover villages than at background points at every buffer radius examined (10–150 km). At 10 km, median cover was 79.2% (IQR 53.7–97.2) among 22 presence rows versus 35.7% (IQR 16.4–59.6) among 10,000 absences. Same-year forest-loss scores were higher at presence sites, with the contrast largest at 10 km and attenuating toward 150 km (Mann–Whitney *p* = 0.14 at 150 km). Fragmentation scores were higher at presence sites at 10, 50, and 150 km.

Log human population within 10 km did not differ detectably (median 9.55 vs 8.73; *p* = 0.33). At 150 km, log population was **lower** at presence sites (median 12.97 vs 14.04; *p* = 0.010). The product of 10 km log population and forest cover (`logpcfc10km`) was markedly higher at presence sites (median 6.42 vs 2.94; *p* = 4.8 × 10⁻⁸).

Precipitation seasonality was lower at presence sites (median 44.8 vs 67.1), as were scaled potential evapotranspiration and temperature seasonality. Elevation and night-time land-surface temperature were similar (*p* > 0.8).

These unadjusted contrasts are consistent with a background of dense forest, recent local forest loss, fragmentation, relatively stable precipitation, and a scale-dependent population pattern. They are **not** BRT partial-dependence curves and should not be read as independent effects after mutual adjustment [16].

**Table 1.** Median (interquartile range) of selected P2-scaled covariates at documented spillover locations versus population-weighted absence locations, 2001–2021.

| Covariate (P2 scaling) | Presence (n) | Median (IQR) | Absence (n) | Median (IQR) | Mann–Whitney *p* |
|---|---:|---|---:|---|---:|
| Forest cover, 10 km (%) | 22 | 79.18 (53.71–97.20) | 10,000 | 35.73 (16.41–59.57) | 4.6 × 10⁻⁷ |
| Forest cover, 100 km (%) | 22 | 83.43 (63.11–89.66) | 10,000 | 35.57 (18.61–56.76) | 7.0 × 10⁻⁷ |
| Forest loss, same year, 10 km | 22 | 28.58 (8.03–73.89) | 10,000 | 2.41 (0.24–12.63) | 5.7 × 10⁻⁵ |
| Forest loss, same year, 100 km | 22 | 18.86 (3.72–46.60) | 10,000 | 6.99 (2.11–18.10) | 0.037 |
| Forest loss, 1 year prior, 10 km | 20 | 22.53 (4.60–98.49) | 9,506 | 2.30 (0.24–11.99) | 8.1 × 10⁻⁵ |
| Forest loss, 2 years prior, 10 km | 18 | 14.54 (4.64–114.13) | 9,000 | 2.39 (0.24–12.10) | 2.4 × 10⁻⁴ |
| Fragmentation, 10 km | 22 | 4.81 (0.00–68.34) | 10,000 | 0.00 (0.00–12.44) | 0.010 |
| Fragmentation, 50 km | 22 | 24.15 (13.65–45.51) | 10,000 | 2.45 (0.00–21.79) | 3.0 × 10⁻⁶ |
| Fragmentation, 150 km | 22 | 20.27 (16.44–41.49) | 10,000 | 7.61 (0.64–22.13) | 1.0 × 10⁻⁴ |
| log(population + 1), 10 km | 22 | 9.55 (6.69–11.38) | 10,000 | 8.73 (7.42–10.01) | 0.33 |
| log(population + 1), 100 km | 22 | 12.92 (10.59–13.86) | 10,000 | 13.41 (12.52–14.45) | 0.061 |
| log(population + 1), 150 km | 22 | 12.97 (11.55–14.50) | 10,000 | 14.04 (13.17–14.95) | 0.010 |
| log population × forest cover, 10 km | 22 | 6.42 (5.49–7.06) | 10,000 | 2.94 (1.39–4.74) | 4.8 × 10⁻⁸ |
| Precipitation seasonality | 22 | 44.80 (33.32–50.75) | 10,000 | 67.10 (51.34–83.02) | 1.8 × 10⁻⁷ |
| PET / 1000 | 22 | 1.37 (1.26–1.54) | 10,000 | 1.63 (1.43–1.85) | 1.6 × 10⁻⁵ |
| Temperature seasonality / 10 | 22 | 58.79 (52.29–77.55) | 10,000 | 97.14 (60.48–130.04) | 8.1 × 10⁻⁵ |
| Elevation (m) | 22 | 485 (445–726) | 10,000 | 560 (370–1,005) | 0.83 |
| Night-time LST | 22 | 293.4 (291.6–293.9) | 10,000 | 293.1 (291.6–294.0) | 0.85 |

Forest-loss scores use the script’s ×10,000 scaling of cover fractions; they are comparable between presence and absence but are not percentages. Tests are unadjusted.

### 3.3 Model-performance figures documented in the scripts and in Telford et al. (2025)

The public repository does **not** contain ROC objects, `model_list` files, or prediction rasters. The following numbers are therefore reported as **documented in `P2 Model fit and CV.R` and/or in the peer-reviewed article**, not as statistics recomputed in this run.

The fitting script hard-codes a confusion table for the multispecies leave-one-year-out predictions:

| | Actual spillover | Actual no spillover |
|---|---:|---:|
| Predicted spillover | 19 | 2,731 |
| Predicted no spillover | 3 | 7,269 |

That table implies sensitivity 19/22 = 86.4% and specificity 7,269/10,000 = 72.7%, which match the published multispecies classification (sensitivity 86.4%, specificity 72.7%) [1]. A nearby comment in the same script listing sensitivity 0.91 is inconsistent with both the table and the paper and should not be used. The script titles an ROC plot “Multi-species Analysis (AUC = 0.88)”; the paper reports leave-one-year-out AUC 0.88 (multispecies) and 0.92 (ZEBOV-only) [1].

A comment block in P2 lists the top 12 predictors from the fitted multispecies ensemble as: log population at 100 km and 10 km; the 10 km population–forest product; fragmentation at 150 km; log population at 150 km; same-year forest loss at 10 km; precipitation seasonality; forest cover at 10 km; temperature seasonality; PET; forest loss two years prior at 10 km; and same-year forest loss at 100 km. The published paper similarly reports human population and forest-related variables among the leading predictors at multiple scales, with opposite signs for population at local versus regional radii, positive associations for forest loss and fragmentation, and higher spillover likelihood where precipitation is more stable [1]. We did not recompute relative-influence medians.

### 3.4 Prediction-grid output that was not reproduced

Script P3 is designed to emit 2021 and 2022 ROR rasters, 1-year change surfaces, and reduced models that drop forest-loss and fragmentation covariates. Those files are absent. Published results that therefore cannot be verified from the public dump include: multispecies 2022 RORs ranging from 0.3 to 32.3, with the top percentile above 10.1 concentrated in DRC (79.1% of top-percentile cells), Cameroon, Gabon, and the Republic of the Congo; ZEBOV-only RORs from 0.2 to 31.2; 1-year ROR increases in 24.9% (multispecies) and 18.4% (ZEBOV-only) of the study area; and the ranking of the 2022 Wangata cell in the top 0.1% of ZEBOV-only RORs (maximum 21.5) versus a modest baseline ROR but a top 0.3% 1-year increase in Mubende [1]. Reduced models that ignored forest change concentrated high ZEBOV-only odds along the equatorial rainforest core, whereas full models raised odds along forest edges and identified West Africa as ROR > 1 only after forest-change covariates were included [1]. Those spatial findings remain attributable to Telford et al. until the grid extracts are shared or independently rebuilt.

---

## 4. Discussion

### 4.1 Comparison with previous ecological studies

Independently tabulated medians in `tdf_clean.csv` align with prior evidence that EVD index cases occur in more forested, more fragmented, and more recently disturbed landscapes than typical background locations. Rulli et al. reported higher fragmentation, but not preferentially higher forest *loss*, around index cases relative to the regional mean [10]. In the Telford table, both same-year loss (especially within 10 km) and fragmentation (especially at 50–150 km) were higher at presence points. That is not a contradiction so much as a difference in contrast: Rulli compared outbreak neighbourhoods with the remainder of West and Central Africa, whereas Telford’s absences are population-weighted random points inside an equatorial, high-rainfall mask [1,10]. Olivero et al. emphasised closed-forest loss within two years along rainforest margins [11]. The deposited series includes same-year, 1-year-prior, and 2-year-prior loss; presence medians remain higher than absence medians at each lag within 10 km, with missing lags for 2001–2002 events.

The scale-dependent population pattern—similar or slightly higher local log population, lower regional log population at presence sites—matches the published BRT partial-dependence description [1] and is compatible with Schmidt et al.’s finding that spillover intensity is not a monotone function of human density [9]. Lower precipitation seasonality at presence sites is consistent with Telford’s report that likelihood was higher where precipitation was more stable [1], and more loosely with work that places Ebola suitability in humid forest climates [8,9].

Judson et al. showed that ecological context differs among Ebolavirus species [4]. The public P2 script fits a multispecies ensemble (Zaire, Sudan, and Bundibugyo together). The ZEBOV-only analysis that, in the paper, placed more weight on forest loss, fragmentation, and low regional population, and produced higher predictions inside dense forest [1], is not implemented as a separate committed script. Species-specific replication would require an explicit second filter on `Strain`.

### 4.2 Implications for surveillance

If annually updated RORs are used operationally, the deposited pipeline is a template rather than a turnkey product. The original authors argued that ranking the extreme tail of RORs or of 1-year increases is more realistic than treating every above-average cell as a site for active search, and that training clinicians in high-ranking patches on recognition and reporting among high-risk occupational groups is a proportionate use of a rare-event model [1]. Analogous climate-informed targeting has been used for Rift Valley fever [20]. Those policy inferences depend on the prediction maps, which this draft does not recreate.

Descriptive results that *can* be checked from the CSV still have a limited operational reading: villages that later reported index cases sat in greener, more fragmented, and more locally deforested neighbourhoods than population-weighted random points. That pattern supports continued attention to forest-edge communities, but it does not identify which 1 km cell will produce the next index case.

### 4.3 Honesty about what this repository does and does not contain

A methods paper that treated hardcoded AUC values as if they had been recomputed, or that described 2022 ROR maps as output of this run, would overstate the public archive. The archive is strongest as (i) a complete covariate table for 24 presence points and 10,000 absences and (ii) an explicit recipe for BRT ensembles and GEE extracts. It is weakest as a bit-reproducible computational artifact: Windows CDC working directories, undeposited rasters, undeposited `.RData` ensembles, retired R spatial packages, and a 100 × leave-one-year-out loop that is expensive to rerun.

Predictive BRT ensembles also should not be interpreted as etiologic models. Coefficients and relative influence can reflect proxies and collinear land-use processes; Westreich and Greenland’s “Table 2 fallacy” warning applies if one reads mutually adjusted landscape terms as independent causal effects [16]. Predictions refer to where an index case is *identified*, which may be the village of residence rather than the exact exposure location [1]. Presence–background designs further assume that absences are locations without spillover rather than locations with undetected spillover; population weighting addresses reporting, not perfect capture of every zoonotic infection.

### 4.4 Limitations of the descriptive re-analysis

The presence sample for 2001–2021 is 22 events. Median contrasts can be estimated; they are sensitive to a handful of West African or densely populated points. Mann–Whitney tests with 10,000 absences are overpowered for tiny shifts and were not multiplicity-adjusted. Absence sampling is already biased toward people; contrasts for population covariates are therefore contrasts against a populated background, not against uniform space. Forest-loss scaling (×10,000) is inherited from the fitting script and is not in natural units. We did not rebuild Hansen or LandScan extracts, so we cannot audit buffer geometry or the 70% fragmentation threshold from the CSV alone.

---

## 5. Conclusion

The public replication repository for Telford et al. (2025) contains a well-documented presence–background table and scripts that specify how annually updated forest, population, and climate layers were turned into a 100-member BRT ensemble and ROR maps for equatorial Africa. Independent analysis of that table shows that documented spillover villages in 2001–2021 had higher local forest cover, higher local forest loss, higher fragmentation, lower precipitation seasonality, and a scale-dependent population signature relative to 10,000 population-weighted background points. Published leave-one-year-out discrimination (multispecies AUC 0.88; sensitivity 86.4%; specificity 72.7%) is internally consistent with the confusion matrix stored in the fitting script but was not recomputed here. Annual relative-odds maps, 1-year change surfaces, and the ranking of the 2022 DRC and Uganda events on those maps remain findings of the original article until the prediction-grid inputs are openly archived or rebuilt. For public-health use, the honest next step is computational reproduction of P2 and P3 from shared rasters, not reuse of undocumented numeric results.

---

## Data and code availability

- Replication repository (upstream): https://github.com/carsontelford/AnnualEbolaPrediction_EID  
- Mirror used for this draft: https://github.com/BoazBaliejukia/AnnualEbolaPrediction_EID  
- Analysis-ready table: `tdf_clean.csv` (24 presence rows, 10,000 absences).  
- Scripts: `P1 Absence point generation.R`, `P2 Model fit and CV.R`, `P3 Prediction grid.R`, plus Google Earth Engine `.txt` extractors for Hansen forest cover, LandScan population, MODIS LST, and FLDAS precipitation seasonality.  
- Peer-reviewed source analysis: Telford et al., *Emerg Infect Dis.* 2025;31(4):689–698. doi:[10.3201/eid3104.241193](https://doi.org/10.3201/eid3104.241193) [1].  
- Preprint of the same analysis: doi:[10.1101/2023.08.29.23294795](https://doi.org/10.1101/2023.08.29.23294795).

Large covariate rasters, fitted `model_list` objects, and 2021–2022 prediction grids are **not** in the public repository.

---

## References

1. Telford CT, Amman BR, Towner JS, Montgomery JM, Lessler J, Shoemaker T. Predictive model for estimating annual Ebolavirus spillover potential. *Emerg Infect Dis.* 2025;31(4):689-698. doi:10.3201/eid3104.241193

2. Rouquet P, Froment JM, Bermejo M, Kilbourn A, Karesh W, Reed P, et al. Wild animal mortality monitoring and human Ebola outbreaks, Gabon and Republic of Congo, 2001–2003. *Emerg Infect Dis.* 2005;11(2):283-290. doi:10.3201/eid1102.040533

3. Leroy EM, Kumulungui B, Pourrut X, Rouquet P, Hassanin A, Yaba P, et al. Fruit bats as reservoirs of Ebola virus. *Nature.* 2005;438(7068):575-576. doi:10.1038/438575a

4. Judson SD, Fischer R, Judson A, Munster VJ. Ecological contexts of index cases and spillover events of different ebolaviruses. *PLoS Pathog.* 2016;12(8):e1005780. doi:10.1371/journal.ppat.1005780

5. Leroy EM, Epelboin A, Mondonge V, Pourrut X, Gonzalez JP, Muyembe-Tamfum JJ, et al. Human Ebola outbreak resulting from direct exposure to fruit bats in Luebo, Democratic Republic of Congo, 2007. *Vector Borne Zoonotic Dis.* 2009;9(6):723-728. doi:10.1089/vbz.2008.0167

6. Centers for Disease Control and Prevention. History of Ebola disease outbreaks. Accessed 2023 as cited by Telford et al. [1]. https://www.cdc.gov/vhf/ebola/history/chronology.html

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

17. Riitters K, Wickham J, O’Neill R, Jones B, Smith E. Global-scale patterns of forest fragmentation. *Conserv Ecol.* 2000;4(2):3.

18. Humanitarian Data Exchange. Health facilities in sub-Saharan Africa. https://data.humdata.org/dataset/health-facilities-in-sub-saharan-africa

19. Golding N, Shearer F. seegSDM: Streamlined functions for species distribution modelling in the SEEG Research Group. https://rdrr.io/github/SEEG-Oxford/seegSDM/man/seegSDM-package.html

20. Anyamba A, Chretien J-P, Small J, Tucker CJ, Formenty PB, Richardson JH, et al. Prediction of a Rift Valley fever outbreak. *Proc Natl Acad Sci U S A.* 2009;106(3):955-959. doi:10.1073/pnas.0806490106
