# Method

Weighted FWL regression on treatment and additive covariates after internal cyclic absorption of cluster (or supplied panel unit) and time FE. Target is a TWFE coefficient, not automatically ATT.

Stata currently uses cluster FE only. TWFE uses native weighted alternating
projections; other models use bounded (12-million estimated-entry) SVD reference
designs. No reghdfe/ftools or cross-language runtime calls. At fixed design the
scalar is linear in Y; deletion/reassignment rebuilds assignment-dependent design.
