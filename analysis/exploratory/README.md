# Exploratory checks

These are the ad-hoc scripts that were run while the analysis was being
developed: one-off cross-checks of a count, a correlation, or a subgroup.
They are kept here for transparency, not because anything in the manuscript
depends on them.

Two things to know before reading them:

- **Nothing reported in the manuscript comes from these scripts.** Every
  number, table and figure in the paper is produced by the numbered scripts
  in the parent directory (`01_` through `11_`), which run in order against
  the tables that the SQL in `../sql/` builds.
- **Several of them will not run as they stand.** They read intermediate
  CSVs by absolute path from the machine the analysis was done on. Those
  files are derived from MIMIC-IV and eICU-CRD and are not redistributable,
  so they are not in this repository and the paths are not portable.

They are included because removing them would curate the record. Read them
as a working notebook, not as part of the pipeline.
