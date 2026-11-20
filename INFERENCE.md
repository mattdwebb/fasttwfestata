# Native inference

OLS cluster CR1, G/(G-1)*(N-1)/(N-K), Student t(G-1). K includes actual additive FE rank and regressors.

## November 2026 support revision

The Stata delete-cluster calculation now recomputes the estimator on surviving support. Disappearance of an adoption cohort by itself does not invalidate a deletion. Every required deletion must yield an estimate for a jackknife standard error.
