# analysis

Python. Each script reads the CSVs exported from the tables that `sql/`
writes, and prints its results while saving every table into `output/`.

Run in order; nothing is shared between them except `config.py`.

| Script | Produces |
|--------|----------|
| `01_main_analysis.py` | primary models in both cohorts, Table 3, the SOFA over-adjustment check, collinearity, discrimination |
| `02_landmark.py` | Table 4, and the arithmetic showing why the primary estimate is conservative |
| `03_cholestatic.py` | alkaline phosphatase, AST/ALT, R ratio, direct bilirubin, in both cohorts |
| `04_sensitivity.py` | eTables 4 to 7: dose–response, competing causes, the definition of HH, the transplant exclusion |
| `05_descriptive_tables.py` | Tables 1, 2 and 5, observation density, survivors versus decedents |
| `06_figures.py` | Figures 1 and 2, TIFF and PNG at 300 dpi |

```bash
pip install -r ../requirements.txt
cd analysis
for f in 0*.py; do python3 "$f"; done
```

Set `HH_DATA` to read the exports from elsewhere and `HH_OUT` to write
somewhere other than `output/`.

## Two conventions worth knowing before reading the code

**`import config as cfg`, never `as C`.** Patsy uses `C()` for categorical
terms inside a formula and evaluates it in the caller's namespace, so binding
`C` to anything else makes every model fail with
`TypeError: 'module' object is not callable`.

**Rounding is half-up, not half-to-even.** `config.rh()` exists because
Python's `round` sends 8.5 to 8, which silently moved a median SOFA and
several LDH values in an earlier version of these tables.
