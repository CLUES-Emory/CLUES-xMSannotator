# Documentation

## Guides

| Document | Description |
|----------|-------------|
| [Pipeline Workflow Reference](xMSannotator_Workflow.md) | Complete 5-stage pipeline reference with scoring formulas, parameter decision guide, and troubleshooting |
| [Input File Formats](xMSannotator_Input_Formats.md) | API-level specification for all input tables, `advanced_annotation()` parameters, and working examples |
| [Input Data Pre-Processing](advanced_annotation_input_formatting.md) | XCMS feature table, sample mapfile, compound database and Step 1 feature group formats with pre-processing steps (feature selection, blank removal, peak table construction) |
| [Stage 5 Output Reference](Stage5_Output_Reference.md) | Column definitions and interpretation guide for Stage5_curated_results.txt |
| [Release notes](../NEWS.md) | Changes by version (1.1.0 onward) |
| [Changelog](CHANGELOG.md) | Dated changes, bug fixes, and new features |

## Example Scripts

| Script | Description |
|--------|-------------|
| [Example Runscript](../inst/scripts/xMSannotator_CLUES_Runscript_Example.R) | Complete single-database annotation workflow with configuration, pre-processing (Step 1 feature groups or fold-change filter), annotation, and result export |
| [Multi-DB Runscript](../inst/scripts/xMSannotator_CLUES_MultiDB_Runscript.R) | SLURM array job R script for multi-database annotation with CLI argument support; reads Step 1 Stage 3 feature groups when present |
| [SLURM Wrapper](../inst/scripts/run_CLUES_xMSannotator.sh) | Bash/SLURM submission script for dual-polarity (pos/neg) batch jobs across multiple compound databases |

## Developer & Legacy Documentation

| Document | Description |
|----------|-------------|
| [Developer Documentation](developer_documentation.md) | Setup instructions, code style (tidyverse + styler), and testing framework |
| [Testing](testing.md) | testthat/patrick testing framework, code coverage tools, and test data |
| [Modifications](modifications.md) | Changes in recetox-xMSannotator vs. [original xMSannotator](https://github.com/kuppal2/xMSannotator) |
| [Refactoring Patterns](refactoring.md) | Design patterns for parameter passing, parallel execution, and data restoration |
| [Possible Issues](possible_issues.md) | Known bugs and issues detected in the original codebase during refactoring |
| [Research Reproducibility](research_reproducibility.md) | Online data dependencies (KEGG, HMDB, ChemSpider) affecting reproducibility |

## Images

- `xMSannotator_Workflow.png` — Pipeline workflow diagram (Figure 1 in README)
- `confidence_levels_v2.png` — Confidence level requirements decision tree (Figure 2 in README)
