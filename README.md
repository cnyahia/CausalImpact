## Causal Impact
This repository includes methods for determining the causal impact of an intervention!

### Description
  * **psm-example.md**: an R markdown file illustrating the use of propensity score matching to infer the causal impact.
  The causal effect under consideration is the impact of labor training on income gain.
  Propensity score matching (psm) controls for confounding variables such as age, race, marital status, etc.
  Precisely, psm creates a balance between the covariate distribution in the treated group and the distribution in the control group.
  Thus, since both groups would have similar covariate distributions, the difference in outcome across groups represents the causal impact (as opposed to impact of confounding variables on outcome).
