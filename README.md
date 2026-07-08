# Individual Tumor Dynamics and Their Link to Survival (Simulated Data)

## Context

This project reproduces, on **fully simulated data**, the methodological approach used
in pharmacometrics to link individual tumor dynamics (longitudinal RECIST-type
measurements) to survival endpoints (PFS/OS) in oncology, and to characterize how drug
exposure translates into tumor response and survival; a central question in clinical
development of anticancer therapies (e.g., Claret et al., Wang et al., *Clin Pharmacol
Ther*, 2009).

No real patient data is used. The goal is to demonstrate methodological understanding,
simulation, individual parameter estimation, joint modeling, and PK/PD-driven tumor
response modeling, as a self-taught learning project, not a replication of any
industry model.

## Notebooks

- **`01_simulation_and_two_stage_model.Rmd`** : simulates 200 patients with individual
  tumor growth/regrowth parameters and survival outcomes, then estimates these
  parameters back from noisy longitudinal data using a classic "two-stage" approach
  (individual nonlinear fits followed by a Cox model). Includes a first attempt with a
  sparse observation design (weak, non-significant signal) and a second, improved design
  (significant signal), illustrating how observation quality affects signal detection.
- **`02_joint_model.Rmd`** : re-analyzes the same (improved-design) data using a joint
  model (`JM` package), estimating the longitudinal trajectory and survival
  simultaneously, and compares the results honestly against the two-stage approach.
- **`03_ADC_PK_PD_extension.Rmd`** : extends the tumor dynamics framework with an
  explicit drug mechanism: a one-compartment PK model and an Emax PD model for a
  simulated antibody-drug conjugate (ADC), linked to a tumor growth inhibition (TGI)
  model and, via a landmark survival analysis, to overall survival. Includes an honest
  discussion of where the model simplifies real ADC pharmacology.

The three notebooks form a coherent progression:

**Notebook 1** → tumor natural history and its link to survival
**Notebook 2** → joint modeling of tumor dynamics and survival
**Notebook 3** → treatment mechanism: ADC PK → PD → TGI → survival

## Models

**Notebooks 1-2** : Tumor growth inhibition (TGI) biexponential model (Wang et al., 2009):

```
SLD(t) = SLD0 * (exp(-g*t) + exp(d*t) - 1)
```

- **g**: tumor regression rate (treatment effect)
- **d**: tumor regrowth rate (resistance)
- **SLD0**: baseline tumor size (Sum of Longest Diameters)

Survival is simulated via a Cox proportional hazards structure in which `g` is
protective and `d` is detrimental.

**Notebook 3** : one-compartment PK model (`C(t) = Dose/V * exp(-kel*t)`), Emax PD
model, and a simplified TGI ODE (`dSLD/dt = g*SLD - Effect(t)*SLD`) driven by plasma ADC
exposure as a surrogate pharmacodynamic driver. Survival is linked through a landmark
analysis using tumor burden at day 60 as a prognostic biomarker, compared against a
no-treatment counterfactual under the same hazard structure.

## Key results

| Approach | Design | cor(g, g_hat) | cor(d, d_hat) | Signal on survival |
|---|---|---|---|---|
| Two-stage | Sparse (6-week scans, 15% noise) | 0.568 | 0.919 | Not significant (p = 0.60) |
| Two-stage | Improved (3-week scans, 8% noise) | 0.830 | 0.990 | Significant (p = 0.007) |
| Joint model | Improved design, simplified linear trend | - | - | Borderline (p = 0.056) |

**Notebooks 1-2 takeaway.** Detecting a true tumor-dynamics/survival association is
highly sensitive to the observation design (scan frequency, measurement precision),
"attenuation by measurement error." The joint model is more statistically rigorous than
the two-stage approach, at the cost of a simplified representation of tumor dynamics, a
genuine methodological trade-off, not a failure of either approach.

**Notebook 3 takeaway.** Individual PK variability (elimination rate) translates into
heterogeneous tumor exposure and response; higher ADC exposure (AUC) is associated with
lower residual tumor burden at day 60, which in turn predicts better post-landmark
survival in this simplified framework.

## Limitations (stated explicitly in each notebook)

- The two-stage approach ignores estimation uncertainty when parameters are reused as
  fixed covariates; the joint model addresses this but uses a simplified linear tumor
  trajectory instead of the full biexponential model.
- The notebook 3 TGI model uses plasma ADC exposure as a *surrogate* driver of
  pharmacodynamic activity, a simplification of the real ADC mechanism (tumor
  distribution, internalization, payload release), and does not model acquired
  resistance or repeated dosing.
- The landmark survival analysis in notebook 3 avoids immortal-time bias but discards
  early events before the landmark and depends on the chosen landmark time.
- Euler integration (1-day step) is used for the TGI ODE; a dedicated ODE solver
  (`deSolve`) would be the more standard implementation for more complex or stiffer
  models.

## Tools

R — packages: `minpack.lm`, `survival`, `survminer`, `nlme`, `JM`, `dplyr`, `ggplot2`

## Author

Juliette Bouli-Mengue 
