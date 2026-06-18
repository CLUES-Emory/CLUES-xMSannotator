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
| 4 | `Module_RTclust` | character | Combined WGCNA module and RT cluster ID (e.g., "46_1" = module 46, RT cluster 1). Features sharing a `Module_RTclust` co-vary in intensity and co-elute. |
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

| Level | Label | Evidence Required |
|:-----:|-------|-------------------|
| 4 | Confirmed | User-verified compound (via `boosted_compounds` parameter) |
| 3 | High | Isotope rows + base adduct(s) + module/RT coherence |
| 2 | Medium | 2+ distinct base adducts + module/RT coherence |
| 1 | Low | Single primary adduct match (M+H or M-H) |
| 0 | None | Insufficient evidence |

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

## See Also

- [Pipeline Workflow Reference](xMSannotator_Workflow.md) — full algorithmic details for each stage
- [Input File Formats](xMSannotator_Input_Formats.md) — input table specifications and `advanced_annotation()` parameters

---

*Document created: 2026-03-14*
*For use with CLUES.xMSannotator v1.0.0*
