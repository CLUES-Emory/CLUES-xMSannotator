# Stage 5 Output Column Reference

> **Note (opt-in output):** As of the 2026 release, `Stage5_curated_results.txt` is produced only when `redundancy_filtering = TRUE`. The canonical final output of the pipeline is now `Stage4b_confidence_levels.txt`, which preserves the full evidence basis behind every confidence value. The column reference below still applies for users who opt in to Stage 5; the columns are identical to Stage 4b plus `MatchCategory`. For confidence-faithful downstream analysis, prefer Stage 4b.

`Stage5_curated_results.txt` contains one row per curated annotation after redundancy filtering — each row is a putative compound identification with a confidence level, chemical score, and match metadata. Rows belonging to losing compounds at multi-match `(mz, time)` features have been deleted.

**Format:** Tab-delimited text file.

---

## Column Reference

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | `compound_id` | character | Compound identifier from the database (e.g., "PEST3133"). If the database provides `compound_id`, those values flow through; otherwise auto-generated as "Formula_1", "Formula_2", etc. |
| 2 | `Confidence` | integer (0–4) | Numeric confidence level from the evidence-based decision tree. See [Interpreting Confidence](#interpreting-confidence). |
| 3 | `score` | numeric | Chemical score from multi-evidence evaluation — incorporates adduct count, correlation strength, isotope confirmation, RT scaling, and confidence-level boost. Higher = stronger evidence. See [Interpreting Score](#interpreting-score). |
| 4 | `Module_RTclust` | character | Combined WGCNA module and RT cluster ID (e.g., "46_1" = module 46, RT cluster 1). Features sharing a `Module_RTclust` co-vary in intensity and co-elute. With `feature_groups` (1.2.0+) the module is the feature group (`as.integer(factor(group))`; for a Step 1 groups file this equals `feature_group`, checked on CLU0127) and the RT cluster is always 1. |
| 5 | `mz` | numeric | Observed mass-to-charge ratio of the detected feature. |
| 6 | `time` | numeric | Retention time in seconds. |
| 7 | `MatchCategory` | character | `"Unique"` = this m/z matched exactly one compound. `"Multiple"` = initially matched multiple compounds; the winner was retained after redundancy filtering. |
| 8 | `theoretical.mz` | numeric | Expected m/z from the compound's monoisotopic mass and adduct formula: `(factor × M + adduct_mass) / charge`. |
| 9 | `delta_ppm` | numeric | Mass error: `10⁶ × |theoretical.mz − mz| / theoretical.mz`, rounded to 2 decimal places. |
| 10 | `Name` | character | Compound name from the database. |
| 11 | `Formula` | character | Molecular formula (e.g., "C41H65NO10"). |
| 12 | `MonoisotopicMass` | numeric | Neutral exact monoisotopic mass in Daltons. |
| 13 | `Adduct` | character | Adduct assignment. For monoisotopic ions: adduct name (e.g., "M+H"). For isotope peaks: adduct with suffix (e.g., "M+ACN+H_[+2]" = second isotope of M+ACN+H). |
| 14 | `mean_int_vec` | numeric | Mean intensity across all samples for this feature. |
| 15 | `MD` | numeric | Mass defect — fractional part of the observed m/z (e.g., 0.5 for m/z 775.4994). |
| 16 | `isotopologue` | character / NA | Specific isotope substitution identified by enviPat (e.g., "18O:1", "13C:2"). NA for monoisotopic peaks or when disabled. |
| 17 | `isotopologue_quality` | character / NA | `"confirmed"` = m/z and intensity match theoretical pattern. `"mz_only"` = m/z matches but intensity does not. NA for monoisotopic peaks. |
| 18 | `Confidence_Level` | character | Human-readable label: "None" (0), "Low" (1), "Medium" (2), "High" (3), "Confirmed" (4). |
| 19 | `feature_id` | character | Custom feature identifier from the input peak table (e.g., "F13164"). Only present if `feature_id_column` was set in `advanced_annotation()`. |

---

## Column Origins

| Origin | Columns |
|--------|---------|
| **Input** (peak table / compound database) | `mz`, `time`, `Name`, `Formula`, `MonoisotopicMass`, `mean_int_vec`, `MD`, `feature_id` |
| **Stage 1** (mass matching) | `compound_id`, `theoretical.mz`, `Adduct` |
| **Stage 1.5** (correlation / clustering) | `Module_RTclust` |
| **Stage 2** (isotope detection) | `Adduct` suffix (e.g., `_[+1]`, `_[+2]`) |
| **Stage 3** (chemical scoring) | `score` |
| **Stage 4** (confidence assignment) | `Confidence`, `Confidence_Level`, `delta_ppm`, `isotopologue`, `isotopologue_quality` |
| **Stage 5** (redundancy filtering) | `MatchCategory` |

---

## Column Name Mapping

Several columns are renamed from their internal names during output formatting (`reformat_annotation_table()` in `R/integration_utils.R`):

| Internal Name | Output Name |
|---------------|-------------|
| `rt` | `time` |
| `expected_mass` | `theoretical.mz` |
| `name` | `Name` |
| `molecular_formula` | `Formula` |
| `monoisotopic_mass` | `MonoisotopicMass` |
| `mean_intensity` | `mean_int_vec` |
| `mass_defect` | `MD` |
| `module` + `rt_cluster` | `Module_RTclust` |
| `cur_chem_score` | `score` |

---

## Interpreting Confidence

| Level | Label | Evidence Required | `schymanski_level` (MultiDB runscript xlsx, 1.6.0) |
|:-----:|-------|-------------------|:-----:|
| 4 | Confirmed | User-verified compound (via `boosted_compounds` parameter) | 1 |
| 3 | High | Isotope rows + base adduct(s) + module/RT coherence | 4a |
| 2 | Medium | 2+ distinct base adducts + module/RT coherence | 4b |
| 1 | Low | Single primary adduct match (M+H or M-H) | 5 |
| 0 | None | Insufficient evidence | 5 |

`schymanski_level` uses the Schymanski et al. 2014 scale (doi:10.1021/es5002105); 4a / 4b split Level 4 (molecular formula) by the evidence used and are not part of the 2014 scheme.

For the complete confidence decision tree, see the [Pipeline Workflow Reference](xMSannotator_Workflow.md#stage-4--confidence-assignment).

## Interpreting Score

| Score Range | Typical Evidence |
|:-----------:|------------------|
| 0 | No corroborating evidence; mass match only |
| 1–10 | Single adduct, low correlation |
| 10–100 | Multiple adducts or weighted adduct with good correlation |
| 100–1,000 | Strong multi-adduct evidence or isotope confirmation |
| 1,000–10,000 | Multiple weighted adducts + isotopes + high correlation |
| >10,000 | Full evidence: weighted adducts + isotopes + confidence boost |

---

## Columns added by the MultiDB runscript (xlsx)

`inst/scripts/xMSannotator_CLUES_MultiDB_Runscript.R` writes `Stage4b_confidence_levels.txt` to an xlsx and adds the columns below. The rows, `Confidence` and `score` are unchanged. Each mode's generated `<mode>-<study>_xMSannotator_README.md` describes them for that run.

| Version | Columns | Joined or computed by |
|---|---|---|
| 1.4.0 | `rt_pred_sec`, `rt_pred_min_sec`, `rt_pred_max_sec`, `rt_pred_source`, `rt_pred_error_sec`, `rt_pred_within_w95`, `p_detected`, `is_rt_standard`, `rt_esi_no_prediction_reason` | `compound_id`, from the reference-library prediction file (empty without one) |
| 1.4.0 | `ms2_status`, `has_ms2`, `ms2_n_peaks`, `ms2_n_files`, `ms2_best_purity`, `ms2_spectrum_method` | `feature_id`, from Step 1 Stage 4 |
| — | database METADATA columns (e.g. `cid`, `pubmed_count`, `patent_count`, `total_references`) | `compound_id`; described in `DATA_DICTIONARY.md` in each database's output folder (1.8.0, copied from the library folder when the bundle has one) |
| 1.5.0 | `isomer_group_n`, `isomer_rank`, `isomer_score`, `isomer_score_noref`, `isomer_rt_term`, `isomer_ion_term`, `isomer_ref_term` | per annotation (`compound_id`, `rt_cluster`), see below; all but `isomer_group_n` blank without a prediction file (1.9.0) |
| 1.6.0 | `schymanski_level` (placed after `Confidence`) | from `Confidence`, see [Interpreting Confidence](#interpreting-confidence) |
| — | feature table QA columns, `feature_group` | `feature_id` |

**Isomer score and rank (1.5.0).** Isomers are annotations with the same `Formula` and the same set of features; xMSannotator gives them identical `Confidence` and `score`. All rows of an annotation (adducts, isotopes) share one value.
- `isomer_rt_term`: -0.5 (RT error / sigma)^2, sigma = window95 / 1.96, both RTs clipped at the void time.
- `isomer_ion_term`: ln(`p_detected`), p limited to [0.001, 0.999].
- `isomer_ref_term`: ln(`total_references` + 1); missing when the count is blank or the row sums several structures.
- Each term is centred on the isomer group's mean; a missing value counts as the mean (0).
- `isomer_score` = 1 x RT + 1 x ionization + 0.5 x references (setting `isomer_score_weights`). `isomer_score_noref` leaves the references out.
- `isomer_rank` is the rank of `isomer_score` in the group (1 = best, ties share the rank). `isomer_group_n` is the group size; with no isomer the score columns are empty and the rank is 1.
- **Without RT / ionization-mode predictions** (`rt_esi_prediction_dir` not set; 1.9.0): the score is not computed. `isomer_rank`, `isomer_score`, `isomer_score_noref` and the three term columns are blank; `isomer_group_n` is still filled. Before 1.9.0 these runs got a rank from the references term alone (RT and ionization terms 0 for every isomer). `run_info.txt` says `isomer_rank: computed` or `not computed`.
- The weights were calibrated on CLU0162 C18 with out-of-fold predictions (truth: Level 1 standards and Level 2a MS2 matches, both well-studied compounds). `total_references` is a popularity prior, not evidence about the feature.

---

## See Also

- [Pipeline Workflow Reference](xMSannotator_Workflow.md) — full algorithmic details for each stage
- [Input File Formats](xMSannotator_Input_Formats.md) — input table specifications and `advanced_annotation()` parameters

---

*Document created: 2026-03-14*
*For use with CLUES.xMSannotator v1.0.0; MultiDB runscript columns and `schymanski_level` added 2026-11-16 (v1.6.0); `DATA_DICTIONARY.md` 2026-11-17 (v1.8.0)*
