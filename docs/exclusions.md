# Cohort construction and exclusions

## Inclusion

An admission enters a cohort when both hold:

1. **Shock** — a diagnosis of septic, cardiogenic or hypovolemic (including
   hemorrhagic and traumatic) shock. See `docs/diagnosis_codes.md`.
2. **Hypoxic hepatitis** — peak alanine aminotransferase >500 U/L at any point
   during the hospital admission.

The ALT threshold follows Raman et al. (Crit Care Med 2002) for severe
ischemic liver injury. Published biochemical definitions of hypoxic hepatitis
range from AST >5× the upper limit of normal (Van den Broecke et al.) to
aminotransferases >20× normal (Henrion et al.), so the manuscript reports the
analysis under four alternative definitions; all give the same conclusion
(`analysis/sensitivity.py`, eTable 6).

## Etiology assignment

For the cross-etiology analysis each admission is assigned to exactly one
etiology by an exclusive hierarchy:

```
septic  >  cardiogenic  >  hypovolemic
```

An admission coded for more than one takes the first that matches. The
validation cohort uses the same hierarchy.

For the secondary hemorrhagic-shock analysis the hierarchy is not applied:
every admission with a hemorrhagic or hypovolemic shock diagnosis is included,
and concurrent other shock types are allowed. The two MIMIC cohorts therefore
differ in size by design, not by discrepancy.

## Exclusion: liver transplantation

Admissions with a liver-transplant procedure (MIMIC-IV) or a diagnosis of
prior or recent liver transplantation (eICU-CRD) are excluded from every
cohort, including the hemorrhagic-shock cohort. The exclusion was specified a
priori, on the grounds that the peri-operative course of transplantation may
itself involve hypotension and marked elevations of aminotransferases and
bilirubin.

### Post-hoc validation of this exclusion

Because an exclusion that removes a tenth of the exposed group can bias an
estimate, the exclusion was checked after the fact
(`sql/07_transplant_exclusion_validation.sql`). Among the 31 excluded
admissions:

- The timing of the AST peak correlated with the date of the transplant
  procedure, Spearman ρ = 0.832.
- In 24 (77.4%) the AST peak fell within two days of the procedure.
- In 5 the AST peak followed the procedure by more than two days.
- In 2 it preceded the procedure by more than two days.

The aminotransferase elevation in these admissions therefore arose from graft
reperfusion rather than from shock. Reinstating the two admissions in which
hypoxic hepatitis appeared to precede transplantation gives an adjusted odds
ratio of 2.38 (95% CI 1.67–3.39) against 2.48 in the primary analysis.

## Competing causes of hyperbilirubinemia

Chronic liver disease, biliary obstruction, hemolysis and rhabdomyolysis are
**not** exclusion criteria. They are identified
(`docs/diagnosis_codes.md`) and handled as sensitivity analyses: patients with
each are excluded in turn, and then all together.

This is a deliberate choice and a limitation. Coded diagnoses under-ascertain
these conditions in both databases, and markedly so in eICU-CRD, where
diagnostic entry is optional — chronic liver disease is recorded in 2.4% of
the validation cohort against 13.0% in MIMIC-IV, and no hemolysis is recorded
at all. The sensitivity analyses therefore demonstrate that the association is
robust to these diagnoses; they do not establish that the conditions were
absent.

## Flow

| | Derivation (MIMIC-IV) | Hemorrhagic (MIMIC-IV) | Validation (eICU-CRD) |
|---|---|---|---|
| Meet shock and ALT criteria | 1,383 | 87 | — |
| Excluded, transplant-related | 31 | 10 | — |
| **Analysed** | **1,352** | **77** | **890** |

The validation cohort is contributed by 112 hospitals, a median of 5 patients
each (IQR 1–9, maximum 52); 49 hospitals contribute at least one patient with
bilirubin ≥10 mg/dL.

## Corrections made during revision

Two errors in an earlier version of this analysis were found and corrected,
and are recorded here because they changed published values:

1. **Windowed laboratory peaks.** Peaks were computed over a restricted window
   rather than the whole admission. Inclusion had always been assessed over the
   whole admission, so the cohort was unaffected, but the reported descriptive
   values were understated. All peaks are now taken over the whole admission.
   This also revealed that peak ALT does not differ between the bilirubin
   groups, which the manuscript now reports.

2. **Category boundaries in the dose–response analysis.** `pd.cut` with its
   default right-closed intervals placed a peak bilirubin of exactly 10.0 in
   the middle category, so the top category held 207 patients rather than the
   209 of the binary exposure. The cut is now half-open on the left
   (`right=False`), and the categorical and binary exposures agree.
