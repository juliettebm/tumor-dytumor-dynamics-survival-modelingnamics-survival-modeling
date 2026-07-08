# Individual Tumor Dynamics and Their Link to Survival (Simulated Data)

## Context

This project reproduces, on **fully simulated data**, the methodological approach used
in pharmacometrics to link individual and multi-lesion tumor dynamics to survival
endpoints (PFS/OS) in oncology, and to characterize how drug exposure translates into
tumor response and survival, a central question in clinical development of anticancer
therapies (e.g., Claret et al., Wang et al., *Clin Pharmacol Ther*, 2009).

No real patient data is used. The goal is to demonstrate methodological understanding,
simulation, individual and multi-lesion parameter estimation, joint modeling, and
PK/PD-driven tumor response modeling, as a self-taught learning project, not a
replication of any industry model.

## Notebooks

- **`01_simulation_and_joint_model.Rmd`** : simulates 200 patients with individual
  tumor growth/regrowth parameters and survival outcomes, then compares two approaches
  to link tumor dynamics to survival: a classic "two-stage" approach (individual
  nonlinear parameter estimation followed by a Cox model) and a joint model (`JM`
  package, simultaneous estimation of the longitudinal trajectory and survival), with
  an honest comparison of both.
- **`02_ADC_PK_PD.Rmd`** : introduces an explicit drug mechanism: a one-compartment PK
  model and an Emax PD model for a simulated antibody-drug conjugate (ADC), single-dose
  administration, linked to a tumor growth inhibition (TGI) model and, via a landmark
  survival analysis, to overall survival.
- **`03_multi_lesion_RECIST_PFS.Rmd`** : extends the framework to multiple individual
  lesions per patient (5-10), with total tumor burden as their sum. Implements explicit
  RECIST 1.1 progression detection (20% relative + 5mm absolute increase from nadir)
  and constructs PFS as a genuine composite endpoint (earliest of progression or
  death), compared against OS — without a drug mechanism.
- **`04_ADC_multi_lesion_PFS_OS.Rmd`** : combines the multi-lesion RECIST framework
  (notebook 3) with the ADC PK/PD mechanism (notebook 2), extended to repeated q3w
  dosing via the pharmacokinetic superposition principle, to compare PFS and OS between
  a treated (ADC) arm and an untreated counterfactual arm.

The notebooks form a coherent progression:

**Notebook 1** → tumor natural history and its link to survival (two-stage vs. joint model)
**Notebook 2** → single-lesion ADC PK/PD/TGI mechanism and landmark survival
**Notebook 3** → multi-lesion RECIST-based PFS vs. OS (no drug mechanism)
**Notebook 4** → multi-lesion RECIST framework combined with repeated-dosing ADC PK/PD,
comparing PFS and OS between treated and untreated arms

## Models

**Notebook 1** : Tumor growth inhibition (TGI) biexponential model (Wang et al., 2009):

```
SLD(t) = SLD0 * (exp(-g*t) + exp(d*t) - 1)
```

- **g**: tumor regression rate (treatment effect)
- **d**: tumor regrowth rate (resistance)
- **SLD0**: baseline tumor size (Sum of Longest Diameters)

Survival is simulated via a Cox proportional hazards structure in which `g` is
protective and `d` is detrimental.

**Notebook 2** : one-compartment PK model (`C(t) = Dose/V * exp(-kel*t)`), Emax PD
model, and a simplified TGI ODE (`dSLD/dt = g*SLD - Effect(t)*SLD`) driven by plasma ADC
exposure as a surrogate pharmacodynamic driver, single-dose administration.

**Notebook 3** : multiple lesions per patient, each following the biexponential TGI
model from notebook 1 with lesion-specific parameter variation; total SLD as their sum;
explicit RECIST 1.1 progression detection.

**Notebook 4** : multiple lesions per patient (as in notebook 3), each driven by a
shared patient-level ADC PK/PD profile (as in notebook 2) extended to repeated q3w
dosing via the superposition principle; RECIST progression and PFS/OS compared between
ADC and no-drug counterfactual arms.

## Key results

| Notebook | Result |
|---|---|
| 1 — Two-stage | cor(g, g_hat) = 0.83, cor(d, d_hat) = 0.99; significant association with survival (p = 0.007) |
| 1 — Joint model | Borderline association (p = 0.056), using a simplified linear tumor trajectory |
| 3 — PFS vs. OS (no drug) | Median PFS (36 weeks) shorter than median OS (52.3 weeks); PFS event rate 98.5% vs. 77% for OS |
| 4 — PFS: ADC vs. No ADC | p < 0.0001, favoring ADC |
| 4 — OS: ADC vs. No ADC | p = 0.021, favoring ADC (smaller effect than on PFS) |

**Notebook 1 takeaway.** The joint model is more statistically rigorous than the
two-stage approach, at the cost of a simplified representation of tumor dynamics, a
genuine methodological trade-off, not a failure of either approach.

**Notebook 3 takeaway.** PFS and OS are genuinely different clocks: progression is
typically detected well before death, so PFS events accumulate faster than OS events in
the same population.

**Notebook 4 takeaway.** A repeated-dosing ADC regimen produces a clear, non-extreme
treatment effect on both PFS and OS, with a more pronounced effect on PFS — consistent
with real oncology trials, where progression is often more sensitive to treatment
effect than overall survival (influenced by subsequent therapy lines and competing
causes of death, not modeled here). An initial single-dose version of the PK model
(as in notebook 2) produced an effect too transient to detect over a long follow-up
window, motivating the switch to repeated dosing.

## Limitations (stated explicitly in each notebook)

- The two-stage approach (notebook 1) ignores estimation uncertainty when parameters
  are reused as fixed covariates; the joint model addresses this but uses a simplified
  linear tumor trajectory instead of the full biexponential model.
- The notebook 2 TGI model uses plasma ADC exposure as a *surrogate* driver of
  pharmacodynamic activity, a simplification of the real ADC mechanism (tumor
  distribution, internalization, payload release), and originally assumed a single
  dose with no re-dosing (addressed in notebook 4).
- The landmark survival analysis in notebook 2 avoids immortal-time bias but discards
  early events before the landmark and depends on the chosen landmark time.
- Lesions (notebooks 3-4) are simulated as varying independently around patient-level
  parameters, a simplification of true inter-lesion heterogeneity; new lesion
  appearance (a RECIST progression criterion) is not modeled, only growth/shrinkage of
  existing lesions.
- In notebook 4, `Emax` was calibrated empirically to produce a realistic, non-extreme
  treatment response range, rather than derived from a mechanistic or literature-based
  value. Dosing continues at fixed intervals regardless of detected progression, which
  would slightly overestimate treatment effect on OS for patients who progress before
  the end of follow-up — discontinuation at progression was left out as a deliberate
  scope decision (would require restructuring the simulation as a sequential loop).
- Euler integration (fixed time step) is used throughout for the TGI ODEs; a dedicated
  ODE solver (`deSolve`) would be the more standard implementation for more complex or
  stiffer models.

## Tools

R - packages: `minpack.lm`, `survival`, `survminer`, `nlme`, `JM`, `dplyr`, `tidyr`, `ggplot2`

## Author

Juliette Bouli-Mengue 
