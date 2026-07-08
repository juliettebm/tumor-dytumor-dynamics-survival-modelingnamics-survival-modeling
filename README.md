# Individual Tumor Dynamics and Their Link to Survival (Simulated Data)

## Context

This project reproduces, on **fully simulated data**, the methodological approach used
in pharmacometrics to link individual tumor dynamics (longitudinal RECIST-type

**Notebooks 1-2 takeaway.** Detecting a true tumor-dynamics/survival association is
highly sensitive to the observation design (scan frequency, measurement precision) —
"attenuation by measurement error." The joint model is more statistically rigorous than
the two-stage approach, at the cost of a simplified representation of tumor dynamics — a
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
  pharmacodynamic activity — a simplification of the real ADC mechanism (tumor
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

