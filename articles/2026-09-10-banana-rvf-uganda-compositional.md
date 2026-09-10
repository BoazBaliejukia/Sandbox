# Banana cultivation and Rift Valley fever outbreaks in southwestern Uganda, 2016–2024: a compositional land-cover analysis from a public grid-cell dataset

**Manuscript type.** Research article draft (unpublished).  
**Provenance.** This draft documents the public analysis release [carsontelford/RVF_Banana_PublicRepo](https://github.com/carsontelford/RVF_Banana_PublicRepo) (Aim 1 of a dissertation project on agricultural land cover and Rift Valley fever in Uganda). Classification performance numbers are taken from the rendered report `Aim1_AnalyzeClassificationbyYearRMD.html` in the companion repository [carsontelford/RCodeDissertation](https://github.com/carsontelford/RCodeDissertation) (knit date 30 October 2025). Descriptive Results were computed directly from the deposited `rvf.csv`. Fitted generalized estimating equation (GEE) objects, substitution odds-ratio tables, inverse-density weights, and manuscript figures are **not** in the public repository. No model coefficients, substitution odds ratios, or effective sample sizes are reported as if they had been estimated in this run.

---

## Abstract

**Background.** Rift Valley fever (RVF) re-emerged in Uganda in 2016 after nearly five decades without confirmed human cases and has since produced recurrent, spatially localized outbreaks, especially in the southwest. Classic East African risk models emphasize anomalous rainfall and flooding. Those meteorological signatures do not fully explain Uganda’s post-2016 pattern, in which events have occurred across seasons and in mixed crop–livestock landscapes. Banana is a dominant perennial crop in southwestern Uganda, but crop-specific land cover has rarely been treated as a compositional exposure.

**Objective.** To describe the study design, data, and intended inference of the public RVF–banana pipeline, and to report only those results that can be obtained from deposited files: the 5 × 5 km grid-cell analysis table and a rendered land-cover classification report.

**Methods.** We reviewed the classification, dataset-assembly, model-fitting, sensitivity, and diagnostic scripts in RVF_Banana_PublicRepo. The analysis is a longitudinal, grid-cell study of southwestern Uganda for 2016–2024. Annual Sentinel-2 composites are classified with year-specific XGBoost models into banana, coffee, ground crops, and non-crop land. Proportions of non-water land enter an isometric log-ratio (ILR) transformation; SuperLearner fits a multivariate inverse-density weight for measured confounders; a binomial GEE with exchangeable working correlation, clustered on grid-cell identity, is specified for outbreak occurrence (`rvfbin`). Substitution contrasts are defined by reallocating a fixed land-cover share from a donor class to a focal class and mapping the change through the ILR coefficients. We tabulated `rvf.csv` and extracted confusion-matrix statistics from the knitted classification HTML. We did not refit SuperLearner or GEE models.

**Results.** The public table contains 27,117 complete cell-years (3,013 cells × 9 years) spanning approximately 29.57–32.49°E and 1.48°S–1.66°N. Ninety cell-years (0.33%) are coded as RVF outbreak-positive, in 84 unique cells (six cells positive in two years). Outbreak cell-years peaked in 2023 (*n* = 36) and 2018 (*n* = 17). Mean banana cover was 22.6% in outbreak cell-years versus 13.7% in non-outbreak cell-years; coffee cover was similar (9.9% vs 8.6%), and ground-crop cover was lower (22.8% vs 27.1%). Outbreak cell-years also had higher mean livestock density and human population density and lower mean CHIRPS annual rainfall. On polygon-held-out test pixels pooled across 2016–2024, collapsed crop versus non-crop classification achieved accuracy 0.9566 (95% CI 0.9563–0.9568), sensitivity 0.9157, specificity 0.9711, and κ = 0.8877; eight-to-ten-class accuracy on the same test pixels was 0.7796. Year-specific collapsed test accuracies ranged from 0.9484 (2022) to 0.9688 (2016). GEE coefficients, weighted-balance diagnostics, and substitution odds ratios are designed outputs of `scripts/03_analysis/Aim1_3bModelFit.R` but are not deposited.

**Conclusions.** The public release is sufficient to describe the cohort, the compositional exposure, and the accuracy of the Sentinel-2 crop classifier. It is not, by itself, a completed etiologic analysis: independent confirmation of whether banana substitutions raise outbreak odds after confounding control still requires a local rerun of the weighted GEE. Unadjusted contrasts in `rvf.csv` are consistent with a banana-enriched, livestock-dense outbreak landscape and should not be read as causal effects.

**Keywords.** Rift Valley fever; banana; compositional data analysis; isometric log-ratio; land cover; Sentinel-2; XGBoost; generalized estimating equations; Uganda.

---

## Outline (IMRaD)

| Section | Content |
|---|---|
| **Title** | Banana cultivation and RVF outbreaks in southwestern Uganda, 2016–2024 |
| **Abstract** | Background, aim, script-based methods, dataset and classification results, missing GEE estimates |
| **Introduction** | RVF ecology; Uganda’s 2016 re-emergence; gap between rainfall-centric models and crop-specific composition; rationale for this pipeline |
| **Methods** | Design; study area and grid; linelist-to-grid outcome; Sentinel-2 XGBoost classification; covariates; ILR, inverse-density weights, GEE, substitutions; what this draft did and did not run |
| **Results** | Structure of `rvf.csv`; outbreak counts by year; unadjusted land-cover and covariate contrasts; classification metrics from the knitted HTML; products that are missing from Git |
| **Discussion** | Comparison with Uganda serosurveys and East African climate models; implications of banana as a perennial exposure; limitations of rare binary outcomes, residual confounding, and an unrun outcome model |
| **Conclusion** | The repository is a reproducible recipe plus a complete analysis table, not a finished set of substitution odds ratios |

---

## 1. Introduction

Rift Valley fever virus (RVFV; family *Phenuiviridae*, genus *Phlebovirus*) is a mosquito-borne zoonosis that causes abortion storms and neonatal death in ruminants and a febrile, occasionally hemorrhagic or encephalitic, illness in people [1]. Human infection is acquired mainly by contact with blood or tissues of infected animals and, less often, by the bite of infected mosquitoes, principally floodwater *Aedes* species that can also maintain the virus transovarially [1,2]. In the Horn of Africa, epizootics have repeatedly followed periods of anomalously heavy rainfall and flooding, a relationship that underpins satellite-based early-warning systems using sea-surface temperature, rainfall, and vegetation indices [3–5].

Uganda’s recent history does not fit that template cleanly. Human RVF was documented in 1968 and then was not laboratory-confirmed again for 48 years [6]. In March 2016, an outbreak in Kabale District, southwestern Uganda, produced four identified human infections and evidence of infection in livestock and in *Aedes* and *Coquillettidia* mosquito pools; whole-genome sequences clustered with the 2006–2007 Kenya-2 clade, consistent with prior introduction rather than a brand-new emergence [6]. A subsequent community serosurvey found IgG seroprevalence as high as 28% in people and 36% in animals in some subcounties, and human seropositivity tracked animal seropositivity, supporting predominantly animal-to-human spillover [7]. Sporadic human and animal outbreaks followed in 2017–2018 and thereafter, concentrated among men with livestock contact and carrying high case-fatality in detected hemorrhagic cases [8,9,9a]. National livestock serosurveys have since estimated overall IgG seroprevalence on the order of 11%, with higher values in cattle than in sheep or goats and with spatial concentrations in southwestern, central, and parts of northeastern Uganda [10,11]. Geostatistical models of the 2017 livestock survey associated seroprevalence with precipitation seasonality, vegetation-index variability, topographic wetness, and human population change [11].

Those studies map *where virus has circulated*. They do not isolate *which agricultural land covers* coincide with the discrete outbreak events that surveillance records after 2016. Land-use change is a plausible modifier of host–vector contact [11,12], and high human density with a large cultivated fraction has been invoked for Kabale specifically [12]. Treating “cropland” as a single class, however, collapses perennial banana stands, coffee, and annual ground crops that differ in canopy, irrigation, and the way livestock are kept among them. Land-cover fractions in a grid cell are compositional: they sum to one, so a naive regression on raw percentages is not identifiable and invites spurious negative correlations [13,14]. Isometric log-ratio (ILR) coordinates map the simplex to Euclidean space and support substitution contrasts (the change in outcome when one part increases at the expense of another) [14].

The public repository analyzed here implements that design for southwestern Uganda, 2016–2024. An accompanying preprint describes the same Aim 1 analysis [15]. The present draft has two narrower aims: (i) to write Methods from the committed scripts, and (ii) to report only Results that can be computed from files that Git actually tracks—the analysis table and a knitted classification report—while stating explicitly what the outcome model is designed to produce but has not deposited.

---

## 2. Methods

### 2.1 Study design and outcome

The intended analysis is a longitudinal, repeated-measures study of 5 × 5 km grid cells. The unit of analysis is the cell-year. The primary outcome in `scripts/03_analysis/Aim1_3bModelFit.R` is `rvfbin`, a binary indicator that at least one geocoded RVF record fell in that cell in that calendar year. A count version (`rvfn`) and a cell-level cumulative count (`cumcases` = sum of `rvfn` over years) are stored in `rvf.csv` but are not the GEE response.

`scripts/02_dataset/Aim1_1PrepFinalDataset.R` reads a raw linelist from `data/raw/rvf_linelist.csv` (not released), spatially joins points to the grid, and aggregates to cell-year counts. Cells with land fraction ≤10% are dropped, and only cells observed in all nine years 2016–2024 are retained. The public `rvf.csv` is the product of that assembly.

This draft treats `rvf.csv` as the analysis-ready table. We did not re-geocode cases or audit the withheld linelist.

### 2.2 Study area and grid

The study polygon is `Aim 1/Shapefiles/SW_UG_StudyArea.shp` (not in the public release). The grid is built at a resolution of 0.0449° (approximately 5 km near the equator) and stored as `grid5km_notclipped.shp`. Coordinates in `rvf.csv` range from 29.569°E to 32.487°E and from 1.485°S to 1.658°N, covering southwestern Uganda, including highland districts around Kabale and adjacent cattle-corridor margins.

### 2.3 Land-cover classification

Annual Sentinel-2 composites [16] are classified separately for each year 2016–2024 (`scripts/01_classification/Aim1_Classification_*LLscript.R`). Training polygons are human-labeled for eight classes: banana, coffee, ground crops, grass, bush, trees, urban, and water. Pixel spectra and focal texture layers (including NDVI and selected bands in dry and wet windows) are extracted with `terra::extract`. Polygons are split 70/30 to train and test so that pixels from the same polygon do not appear in both sets. Class imbalance is addressed with inverse-frequency weights. `caret::train` fits `xgbTree` models [17] with objective `multi:softprob`, a custom collapsed-crop F1(β = 0.5) metric, and a small grid over `nrounds` ∈ {1000, 1200, 1400}, `eta` ∈ {0.10, 0.12, 0.14}, `gamma` ∈ {2, 3}, `max_depth` = 5, `colsample_bytree` = 0.6, `min_child_weight` = 3, and `subsample` = 0.5.

Class-probability mosaics are temporally smoothed (C++ helpers in `scripts/04_sensitivity/my_functions.cpp`). Water is set to missing before proportions are computed, so banana, coffee, ground-crop, and non-crop fractions are shares of *non-water land*. `exactextractr` sums binary class rasters to the 5 km grid (`scripts/02_dataset/Aim1_1PrepFinalDataset.R`). Variables `bananapct_marg`, `coffeepct_marg`, `gcpct_marg`, and `noncroppct` are the four-part composition used in the main ILR model; they sum to one in the deposited table.

### 2.4 Environmental and population covariates

Covariates were extracted in Google Earth Engine and summarized to the grid as cell means (or sums where noted):

| Domain | Script path / product name | Role in the outcome model |
|---|---|---|
| Precipitation | `GEE_CovarRasters/Rain_chirps` (CHIRPS annual sum, coefficient of variation, cell-specific z-score anomaly) [18] | Confounder |
| Alternative rain / temperature | FLDAS rain and temperature stacks | Available in the table; CHIRPS terms enter the weight model |
| Topography | Elevation (DEM), slope in degrees, topographic wetness index \( \log[(a+1)/\tan\beta] \) | Confounder |
| Soils | OpenLandMap clay (and sand, haplic planosols) | Clay enters the weight model |
| Livestock | FAO Gridded Livestock of the World area-weighted cattle, sheep, and goats (`livestock_aw`) | Confounder (log) |
| Human population | Gridded Population of the World density (`gpwPD`); LandScan, WorldPop, GHS-POP, and Meta layers are also stored | `gpwPD` (log) enters the weight model |
| Access | Distance to highway, water, and parks | Distance-to-park is censored at 0 if negative; not in the default SuperLearner set |

The SuperLearner nuisance model in `Aim1_3bModelFit.R` uses year, `slopeDegrees`, `livestock_aw_log`, `gpwPD_log`, `dem_log`, `clay`, `twi`, `chirpssum_log`, `chirpsCoV`, and `chirpssum_anomz`.

### 2.5 Analytic sample restrictions coded in the model script

Before weighting, the script (i) applies a Rogan–Gladen-style correction [19] to banana, coffee, and ground-crop proportions using year-specific sensitivity and specificity constants hard-coded in the script; (ii) drops cells with `slopeDegrees` ≥ 28; and (iii) again requires nine observed years per cell. The hard-coded classification constants are listed in Table 3. They are *inputs* to a sensitivity-adjusted composition, not outputs of this draft’s computations.

### 2.6 Compositional exposure, inverse-density weights, and GEE

Four-part compositions \( (p_{\text{banana}}, p_{\text{coffee}}, p_{\text{gc}}, p_{\text{noncrop}}) \) are closed and mapped to three ILR coordinates with `compositions::ilr` [14]. A SuperLearner ensemble of `SL.glm` and `SL.gam` [20], fit in parallel with seed 2343, predicts each ILR coordinate from the confounder set. Residuals are treated as multivariate normal; the inverse-density weight is the ratio of the marginal ILR density to the conditional residual density (`mvtnorm::dmvnorm`). Weights are truncated at the 1st and 98th percentiles. Weighted correlations and quintile standardized mean differences are computed as balance checks.

The outcome model is

\[
\operatorname{logit}\{P(\texttt{rvfbin}_{it}=1)\} = \beta_0 + \beta_1\,\mathrm{ILR}1_{it} + \beta_2\,\mathrm{ILR}2_{it} + \beta_3\,\mathrm{ILR}3_{it},
\]

fit with `geepack::geeglm`, binomial family, exchangeable working correlation, cluster `id = ID`, and weights `sw_joint_trunc` [21]. Substitution log odds ratios are linear contrasts \( \Delta z^{\top}\beta \) with delta-method standard errors from the GEE robust covariance, for 5- and 10-percentage-point reallocations from each donor class (or a proportional split across the other three) to each focal class. A simpler banana-only GEE, weighted by a univariate inverse-density weight on `bananapct_marg`, is specified as a secondary “take-home” model.

### 2.7 Sensitivity analysis

`scripts/04_sensitivity/Aim1_SensitivityAnalysis.R` repeatedly samples class labels from the eight-class probability stacks, reapplies temporal smoothing, rebuilds grid proportions, and refits the weighted ILR-GEE. Those iterations require mosaics and a `dataset_withoutCrops.csv` that are not in the public tree.

### 2.8 What this draft computed

We (i) profiled `rvf.csv` with the Python standard library (row counts, outbreak tallies, means, standard deviations, and percentiles); (ii) extracted printed `caret::confusionMatrix` blocks from the knitted HTML in RCodeDissertation; and (iii) recorded script constants. We did not install SuperLearner, geepack, or compositions, and we did not produce odds ratios.

---

## 3. Results

### 3.1 Structure of the public analysis table

`rvf.csv` has 27,117 rows and 88 analytic fields (plus a row index): 3,013 cells observed in each year 2016–2024. Banana, coffee, ground-crop, and non-crop proportions sum to 1.00 on every row. Mean land fraction is 0.918 (minimum 0.101, matching the >10% land filter). Nine cell-years have slope ≥28°; dropping them and requiring a complete nine-year panel leaves 27,108 rows and the same 90 outbreak-positive cell-years.

### 3.2 Outbreak occurrence

Ninety cell-years (0.33%) have `rvfbin = 1`, in 84 unique cells. Six cells are positive in two different years (maximum repeats = 2). Summing each cell’s `cumcases` once yields 142 case-records assigned to those 84 cells (maximum 12 in a single cell). Outbreaks occur in every study year (Table 1), with marked peaks in 2018 and 2023.

**Table 1.** RVF-positive 5 km cell-years in `rvf.csv`.

| Year | Cells | Outbreak-positive cell-years | Percent |
|---:|---:|---:|---:|
| 2016 | 3,013 | 3 | 0.10 |
| 2017 | 3,013 | 4 | 0.13 |
| 2018 | 3,013 | 17 | 0.56 |
| 2019 | 3,013 | 3 | 0.10 |
| 2020 | 3,013 | 5 | 0.17 |
| 2021 | 3,013 | 4 | 0.13 |
| 2022 | 3,013 | 4 | 0.13 |
| 2023 | 3,013 | 36 | 1.19 |
| 2024 | 3,013 | 14 | 0.46 |
| **Total** | **27,117** | **90** | **0.33** |

### 3.3 Unadjusted land cover and covariates

Table 2 summarizes the four-part composition and the confounders that enter the weight model. These are crude contrasts. They are not substitution effects and are not adjusted for spatial clustering or confounding.

**Table 2.** Means (SD) in `rvf.csv` by outbreak status. Land-cover terms are proportions of non-water land.

| Variable | Non-outbreak cell-years (*n* = 27,027) | Outbreak cell-years (*n* = 90) |
|---|---:|---:|
| Banana proportion | 0.137 (0.146) | 0.226 (0.192) |
| Coffee proportion | 0.086 (0.097) | 0.099 (0.108) |
| Ground-crop proportion | 0.271 (0.216) | 0.228 (0.211) |
| Non-crop proportion | 0.506 (0.322) | 0.447 (0.301) |
| All crop proportion | 0.494 (0.322) | 0.553 (0.301) |
| CHIRPS annual rainfall (mm) | 1,170 (187) | 1,069 (150) |
| CHIRPS coefficient of variation | 1.694 (0.222) | 1.905 (0.281) |
| CHIRPS rainfall z-anomaly | −0.000 (0.943) | 0.084 (0.823) |
| GPW human population density | 195 (169) | 279 (272) |
| FAO livestock density (area-weighted) | 4,110 (3,046) | 5,370 (2,735) |
| Elevation (m) | 1,300 (293) | 1,411 (263) |
| Slope (degrees) | 5.31 (4.20) | 6.01 (3.77) |
| Topographic wetness index | 6.84 (3.40) | 6.57 (2.82) |
| Clay (%) | 37.95 (3.22) | 37.78 (3.08) |

Banana cover is higher, and ground-crop and non-crop cover lower, where outbreaks were recorded. Coffee differs little. Outbreak cell-years are wetter in *relative intra-annual variability* (higher CHIRPS CoV) but *drier in annual total*, and they sit in denser livestock and human populations at slightly higher elevation. Median banana cover is 0.149 in outbreak cell-years versus 0.086 otherwise; the banana contrast is therefore not driven only by a few extreme cells.

### 3.4 Classification performance (knitted report, 30 October 2025)

The HTML report evaluates XGBoost predictions on extracted training-polygon pixels, not on the 5 km analysis grid. Train and test pixels are split by polygon. Collapsed “crop versus non-crop” uses the script’s binary `pred_probxgb`. Eight-to-ten-class tables use `pred_classxgb` (class0–class9 labels in the report; the classification scripts themselves train eight named classes).

**Pooled years, collapsed crop vs non-crop**

| Split | Accuracy (95% CI) | Sensitivity | Specificity | κ | Balanced accuracy |
|---|---:|---:|---:|---:|---:|
| Train | 0.9934 (0.9933–0.9934) | 0.9865 | 0.9956 | 0.9822 | 0.9911 |
| Test | 0.9566 (0.9563–0.9568) | 0.9157 | 0.9711 | 0.8877 | 0.9434 |

The pooled test confusion matrix is TP 853,900; FP 75,680; FN 78,603; TN 2,543,023. Pooled ten-class accuracy is 0.9754 on the training pixels and **0.7796** (95% CI 0.7791–0.7800) on the test pixels—the drop is the relevant generalization figure for species-level labels.

**Year-specific collapsed test accuracy** (crop vs non-crop): 2016 0.9688; 2017 0.9563; 2018 0.9522; 2019 0.9614; 2020 0.9552; 2021 0.9548; 2022 0.9484; 2023 0.9596; 2024 0.9519. Corresponding test sensitivities range from 0.8902 (2022) to 0.9343 (2020); specificities range from 0.9627 (2020) to 0.9824 (2016).

**Table 3.** Year-specific sensitivity and specificity hard-coded in `Aim1_3bModelFit.R` for Rogan–Gladen adjustment of banana, coffee, and ground-crop proportions. These constants are documented in the model script; they are not recomputed here.

| Year | Banana Se / Sp | Coffee Se / Sp | Ground crops Se / Sp |
|---:|---|---|---|
| 2016 | 0.90 / 0.98 | 0.84 / 0.98 | 0.70 / 0.96 |
| 2017 | 0.88 / 0.98 | 0.85 / 0.98 | 0.68 / 0.96 |
| 2018 | 0.92 / 0.98 | 0.84 / 0.99 | 0.79 / 0.97 |
| 2019 | 0.94 / 0.98 | 0.86 / 0.99 | 0.78 / 0.97 |
| 2020 | 0.95 / 0.99 | 0.86 / 0.98 | 0.82 / 0.97 |
| 2021 | 0.95 / 0.99 | 0.87 / 0.99 | 0.79 / 0.98 |
| 2022 | 0.92 / 0.99 | 0.88 / 0.99 | 0.81 / 0.97 |
| 2023 | 0.93 / 0.99 | 0.90 / 0.99 | 0.85 / 0.98 |
| 2024 | 0.93 / 0.99 | 0.90 / 0.99 | 0.82 / 0.97 |

Ground-crop sensitivity is systematically lower than banana or coffee, which is consistent with greater spectral confusion among annual crops, grass, and bush. Stochastic (`*_stoch`) land-cover columns are present in `rvf.csv` for a single stored draw; the full Monte Carlo loop is not.

### 3.5 Intended model products that are not in Git

`Aim1_3bModelFit.R` is written to print GEE summaries, write `Aim 1/Model fit results/MagnitudeBaselineOR_curvesDf.csv`, and save figures under `Aim 1/ManuscriptFiles/` (confounder scatter grids, boxplots, positivity heatmaps, pre/post-weighting correlation and SMD plots, and substitution odds-ratio curves for banana, coffee, ground crops, and non-crop land). None of those paths exist in the public repository. The script’s later banana-only block also references objects (`fitglmmTMB`, `bananaModContainer`) that are not created in the committed file, so that appendix is not a runnable standalone analysis as deposited.

We therefore do **not** report substitution odds ratios, ILR coefficients, weighted correlations, SMDs, or effective sample size.

---

## 4. Discussion

### 4.1 What the deposited table can and cannot support

The public grid is a complete nine-year panel with a rare binary outcome (90 events). Unadjusted, outbreak cell-years are banana-richer, slightly coffee-similar, ground-crop-poorer, livestock-denser, more populous, higher, and drier in annual CHIRPS total than non-outbreak cell-years. That pattern is compatible with the repository’s scientific hypothesis—that banana stands, not cropland in general, mark outbreak landscapes—but it is also compatible with confounding by livestock, highland settlement, and surveillance intensity. The ILR-GEE with inverse-density weights is the authors’ answer to that confounding [14,20,21]. Because the fitted objects are absent, this draft cannot claim that the banana contrast survives adjustment.

A June 2026 preprint of the same Aim 1 analysis reports that, after weighting, banana was the only class whose substitutions increased outbreak odds, including a 10-percentage-point shift into banana at low baseline banana cover (odds ratio 1.64, 95% CI 1.17–2.29) and a banana-only model odds ratio of 1.21 (95% CI 1.02–1.43) per 10-percentage-point increase [15]. Those numbers are literature, not results of this run. They should be treated as unreproduced until `Aim1_3bModelFit.R` is executed against `rvf.csv` on a machine with SuperLearner, geepack, and compositions.

### 4.2 Comparison with existing RVF epidemiology in Uganda and East Africa

The 2016 Kabale investigation and follow-up serosurveys established that RVFV was already circulating in southwestern livestock and people before—and after—the first recognized human cluster [6,7,11]. Tumusiime and colleagues mapped national livestock IgG risk to lower precipitation seasonality, haplic planosols, and, contrary to a simple “more cattle, more virus” story, *lower* cattle density in some specifications [10]. Arinaitwe and colleagues, reviewing livestock outbreak records for 2013–2022, found seropositivity concentrated in the central and southwestern cattle corridor and noted land-use pressure around Kabale [12]. The present grid sits inside that southwestern focus. The crude CHIRPS contrast—lower annual totals but higher intra-annual variability in outbreak cell-years—differs from the canonical Kenyan flood signal [3–5] and is closer to the preprint’s claim that Uganda’s post-2016 events are poorly aligned with meteorological anomalies [15]. That does not refute a rainfall mechanism; CHIRPS cell means can miss local flooding, and a single annual sum is a coarse proxy for the short rains that hatch *Aedes* eggs [2,18].

Human infection in the Ugandan case series remains overwhelmingly associated with slaughter, butchery, and herding [7–9,9a]. A land-cover association, if confirmed, would operate upstream of that contact: banana landscapes may concentrate people, smallholder livestock, and residual moisture without being a direct transmission route. The scripts do not measure mosquito abundance, farm practices, or market chains, so any mechanism remains speculative.

### 4.3 Why composition matters

Cropland as a single percentage cannot distinguish a banana substitution from a ground-crop substitution. ILR plus g-computation-style reallocations are the appropriate language for “what happens if 10% of a cell that is now coffee becomes banana” [13,14]. The cost is interpretability: ILR coefficients are not themselves “banana effects.” The repository correctly translates them back to pairwise substitutions. Until those translations are rerun, Table 2 should not be paraphrased as “banana causes outbreaks.”

### 4.4 Classification error

Collapsed crop versus non-crop performance on held-out polygons is strong (test accuracy 0.957, κ = 0.89). Species-level accuracy is materially weaker (test accuracy 0.780). Ground-crop sensitivity in the hard-coded adjustment table is the weak point (0.68–0.85). Misclassifying banana as ground crop, or the reverse, would bias the very substitutions the GEE is meant to compare. The stochastic pipeline is the authors’ response; only one stochastic draw is stored in `rvf.csv`. Pixel-level accuracy also overstates grid-level accuracy if errors are spatially correlated inside a 5 km cell.

### 4.5 Limitations

Several limitations are structural. First, `rvfbin` records *detected* events. Viral hemorrhagic fever surveillance expanded after 2010 [6,9]; spatial variation in testing can mimic spatial variation in risk. Second, 90 events limit how many confounders can be modeled in the outcome equation; the design puts those confounders in the weight model instead, which is reasonable but sensitive to positivity and to SuperLearner specification (only GLM and GAM are in the committed library, despite ranger and xgboost being loaded). Third, livestock and GPW layers are not truly annual for every year; using a static or two-date mean as a time-varying confounder can leave residual confounding. Fourth, slope filtering and the 10% land rule change the target population: lakeshores and very steep highland faces are under-represented. Fifth, the banana-only appendix in the public model script is internally inconsistent and should not be executed as committed. Sixth, this draft did not refit the outcome model, so it cannot confirm or refute the preprint’s odds ratios [15].

---

## 5. Conclusion

RVF_Banana_PublicRepo is a handoff-ready specification for a compositional land-use analysis of RVF in southwestern Uganda: a complete 2016–2024 grid, a documented Sentinel-2 XGBoost classifier with strong collapsed-crop accuracy, and a fully written ILR–inverse-density-weight–GEE pipeline. The deposited table shows a rare, recurrent outbreak process whose crude geography is banana-enriched and livestock-dense. That is a description, not an effect estimate. A paper that claims banana substitutions raise outbreak odds still depends on reproducing the weighted GEE and the sensitivity loop that the public tree points to but does not contain.

---

## References

1. Pepin M, Bouloy M, Bird BH, Kemp A, Paweska J. Rift Valley fever virus (*Bunyaviridae: Phlebovirus*): an update on pathogenesis, molecular epidemiology, vectors, diagnostics and prevention. *Vet Res*. 2010;41(6):61. doi:10.1051/vetres/2010033

2. Linthicum KJ, Anyamba A, Tucker CJ, Kelley PW, Myers MF, Peters CJ. Climate and satellite indicators to forecast Rift Valley fever epidemics in Kenya. *Science*. 1999;285(5426):397-400. doi:10.1126/science.285.5426.397

3. Anyamba A, Chretien JP, Small J, Tucker CJ, Formenty P, Richardson JH, et al. Prediction of a Rift Valley fever outbreak. *Proc Natl Acad Sci USA*. 2009;106(3):955-959. doi:10.1073/pnas.0806490106

4. Anyamba A, Linthicum KJ, Small J, Britch SC, Pak E, de La Rocque S, et al. Prediction, assessment of the Rift Valley fever activity in East and Southern Africa 2006–2008 and possible vector control strategies. *Am J Trop Med Hyg*. 2010;83(2 Suppl):43-51. doi:10.4269/ajtmh.2010.09-0289

5. Linthicum KJ, Bailey CL, Davies FG, Tucker CJ. Detection of Rift Valley fever viral activity in Kenya by satellite remote sensing imagery. *Science*. 1987;235(4796):1656-1659. doi:10.1126/science.3823909

6. Shoemaker TR, Nyakarahuka L, Balinandi S, Ojwang J, Tumusiime A, Mulei S, et al. First laboratory-confirmed outbreak of human and animal Rift Valley fever virus in Uganda in 48 years. *Am J Trop Med Hyg*. 2019;100(3):659-671. doi:10.4269/ajtmh.18-0732

7. Nyakarahuka L, de St. Maurice A, Purpura L, Ervin E, Balinandi S, Tumusiime A, et al. Prevalence and risk factors of Rift Valley fever in humans and animals from Kabale district in Southwestern Uganda, 2016. *PLoS Negl Trop Dis*. 2018;12(5):e0006412. doi:10.1371/journal.pntd.0006412

8. Birungi D, Aceng FL, Bulage L, Nkonwa IH, Mirembe BB, Biribawa C, et al. Sporadic Rift Valley fever outbreaks in humans and animals in Uganda, October 2017–January 2018. *J Environ Public Health*. 2021;2021:8881191. doi:10.1155/2021/8881191

9. Nyakarahuka L, Balinandi S, Mulei S, Kyondo J, Tumusiime A, Klena J, Lutwama JJ, Shoemaker TR. Ten outbreaks of Rift Valley fever in Uganda 2016–2018: epidemiological and laboratory findings. *Int J Infect Dis*. 2019;79(Suppl 1):4. doi:10.1016/j.ijid.2018.11.029

9a. Nyakarahuka L, Whitmer S, Klena J, Balinandi S, Talundzic E, Tumusiime A, et al. Detection of sporadic outbreaks of Rift Valley fever in Uganda through the national viral hemorrhagic fever surveillance system, 2017–2020. *Am J Trop Med Hyg*. 2023;108(5):995-1002. doi:10.4269/ajtmh.22-0410

10. Tumusiime D, Isingoma E, Tashoroora OB, Ndumu DB, Bahati M, Nantima N, et al. Mapping the risk of Rift Valley fever in Uganda using national seroprevalence data from cattle, sheep and goats. *PLoS Negl Trop Dis*. 2023;17(5):e0010482. doi:10.1371/journal.pntd.0010482

11. Telford C, Nyakarahuka L, Waller L, Kitron U, Shoemaker T. Geostatistical modeling and prediction of Rift Valley fever seroprevalence among livestock in Uganda. *Am J Trop Med Hyg*. 2023;108(4):712-721. doi:10.4269/ajtmh.22-0555

12. Arinaitwe E, Atuhaire DK, Hasahya E, Nakanjako GK, Mwebe R, Nizeyimana G, et al. Spatial and temporal analysis of Rift Valley fever outbreaks in livestock in Uganda: a retrospective study from 2013 to 2022. *BMC Vet Res*. 2025;21:375. doi:10.1186/s12917-025-04825-6

13. Aitchison J. The statistical analysis of compositional data. *J R Stat Soc Series B*. 1982;44(2):139-160. doi:10.1111/j.2517-6161.1982.tb01195.x

14. Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barceló-Vidal C. Isometric logratio transformations for compositional data analysis. *Math Geol*. 2003;35(3):279-300. doi:10.1023/A:1023818214614

15. Telford C, Nyakarahuka L, Baluku J, Mutesi J, Song C, Boyce R, Emch M, Edwards J, Shoemaker T, Lessler J. Crop-specific agricultural land use and Rift Valley fever outbreaks in Uganda: a longitudinal compositional analysis. *bioRxiv*. 2026. doi:10.64898/2026.06.24.734419

16. Drusch M, Del Bello U, Carlier S, Colin O, Fernandez V, Gascon F, et al. Sentinel-2: ESA’s optical high-resolution mission for GMES operational services. *Remote Sens Environ*. 2012;120:25-36. doi:10.1016/j.rse.2011.11.026

17. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. In: *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*. New York: ACM; 2016. p. 785-794. doi:10.1145/2939672.2939785

18. Funk C, Peterson P, Landsfeld M, Pedreros D, Verdin J, Shukla S, et al. The climate hazards infrared precipitation with stations—a new environmental record for monitoring extremes. *Sci Data*. 2015;2:150066. doi:10.1038/sdata.2015.66

19. Rogan WJ, Gladen B. Estimating prevalence from the results of a screening test. *Am J Epidemiol*. 1978;107(1):71-76. doi:10.1093/aje/107.1.71

20. van der Laan MJ, Polley EC, Hubbard AE. Super Learner. *Stat Appl Genet Mol Biol*. 2007;6:Article 25. doi:10.2202/1544-6115.1309

21. Liang KY, Zeger SL. Longitudinal data analysis using generalized linear models. *Biometrika*. 1986;73(1):13-22. doi:10.1093/biomet/73.1.13

---

## Data and code

- Public dataset and scripts: https://github.com/carsontelford/RVF_Banana_PublicRepo (`rvf.csv`; `scripts/01_classification` through `scripts/05_diagnostics`).
- Companion classification report: https://github.com/carsontelford/RCodeDissertation (`Aim1_AnalyzeClassificationbyYearRMD.html`, 30 October 2025).
- Intended manuscript citation specified in the public README: Telford C, Nyakarahuka L, Baluku J, Mutesi J, Song C, Boyce R, Emch M, Edwards J, Shoemaker T, Lessler J. Banana cultivation and Rift Valley fever outbreaks in Uganda: a longitudinal compositional land-use analysis.
- Raw Sentinel-2 stacks, training polygons, shapefiles, the RVF linelist, Earth Engine rasters, fitted SuperLearner/GEE objects, and substitution-curve CSVs are not in the public release.

**Draft note.** Descriptive statistics in Tables 1–2 were computed from `rvf.csv` on 10 September 2026. Classification statistics in Section 3.4 were copied from the knitted HTML, not recomputed from withheld pixel CSVs. No GEE was fitted for this draft.
