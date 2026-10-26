# SDM_inputprep

This workflow reconstructs the analysis-ready input file from occurrence and climate data obtained directly by each authorized user. This repository contains neither National Forest Inventory records nor Agro-Meteorological Grid Square Data (AMGSDS).

#### Data that users must obtain

1. Obtain the tree-record and plot files from the Forestry Agency of Japan’s National Forest Inventory under the applicable terms of use.
2. Register to use the Agro-Meteorological Grid Square Data System (AMGSDS) provided by the National Agriculture and Food Research Organization (NARO).
3. Obtain the current IDCS/MFA-compatible `AMD_Tools4.py` from the authenticated AMGSDS user documentation and place it in this directory. Do not commit or redistribute this file unless its terms of use explicitly permit redistribution.

Do not store user IDs, passwords, access tokens, device codes, or MFA information in this repository.

#### Build plot-level occurrence data

```bash
Rscript Build_occurrence_input.R \
  /path/to/06_様式３－大円.csv \
  /path/to/01_プロット.csv \
  data_work
```

The plot file defines the master sampling frame. Consequently, forest plots in which neither focal species was recorded are retained as valid double absences. Japanese species names are matched exactly after whitespace normalization:

- `ミズナラ` → `Qmon`
- `コナラ` → `Qser`

This step generates:

- **`Qmon_Qser_occurrence.csv`**: five-column input file used for climate retrieval
- **`occurrence_audit.csv`**: audit of record counts, file encodings, coordinates, and data joins

#### Retrieve and aggregate climate data

```bash
python Fetch_amgsds_climate.py \
  --input data_work/Qmon_Qser_occurrence.csv \
  --outdir data_work \
  --start-year 2016 \
  --end-year 2025 \
  --min-coverage 0.95 \
  --min-valid-years 10
```

The annual climate file is updated after each plot has been retrieved successfully. If the process stops because of a network, session, or authentication error, rerun the same command. Plots already completed will be detected and skipped.

Use `--no-resume` only when intentionally rebuilding the annual climate file from the beginning.

#### Climate aggregation

All seven candidate climate variables are retained in the intermediate files and in `Quercus_climate.csv`. This allows collinearity screening and alternative variable selection to be reproduced downstream.

| Output | AMGSDS daily element | Annual aggregation | 2016–2025 aggregation |
|---|---|---|---|
| `AMT` | `TMP_mea` | Mean | Mean of annual values |
| `Tmax` | `TMP_max` | Maximum | Mean of annual values |
| `Tmin` | `TMP_min` | Minimum | Mean of annual values |
| `APCP` | `APCP` | Sum | Mean of annual values |
| `SD` | `SD` | Maximum | Mean of annual values |
| `SSD` | `SSD` | Sum | Mean of annual values |
| `RH` | `RH` | Mean | Mean of annual values |

An annual value is accepted only when at least 95% of its expected daily observations are finite. With the default setting `--min-valid-years 10`, a 2016–2025 climate value is accepted only when all ten annual values are valid.

Accordingly, `Tmax` and `Tmin` represent the means of the annual extremes rather than the single maximum and minimum values recorded across the entire decade.

#### Final analysis input

The workflow produces `Quercus_climate.csv`, containing one row per retained forest plot and the following 12 columns:

```text
plot_id,Qmon,Qser,longitude,latitude,AMT,Tmax,Tmin,APCP,SD,SSD,RH
```

Plots with invalid coordinates or missing values in any of the seven climate variables are excluded from this file and recorded separately in `excluded_plots.csv`.

#### Generated files

- **`annual_climate_2016_2025.csv`**: plot × year climate values and daily-data coverage
- **`climate_2016_2025.csv`**: 2016–2025 summaries of all seven candidate variables
- **`Quercus_climate.csv`**: 12-column input file for the SDM analyses
- **`excluded_plots.csv`**: excluded plots and their exclusion reasons
- **`Quercus_climate_audit.csv`**: final row counts and species-occurrence counts
- **`aggregation_settings.json`**: machine-readable climate aggregation settings

The `.gitignore` excludes source data, generated derivatives, `AMD_Tools4.py`, local authentication files, Python virtual environments, and cache files.

Users are responsible for reviewing and complying with the current terms of use of both data providers. Both datasets should be cited appropriately in publications based on the reconstructed input.