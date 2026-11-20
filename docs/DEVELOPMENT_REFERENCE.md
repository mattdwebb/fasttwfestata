> Historical development notes retained for provenance. For current behavior, installation and examples, use the [main user guide](../README.md). Old audit statements below may describe superseded behavior.

# fasttwfestata

Native Stata/Mata static TWFE, version 0.1.0.
PARTIAL: frozen-historical comparisons pass on StataNow 19.5; not full feature
parity or fresh official estimator-software certification.

Install with `net install fasttwfe, from("LOCAL_REPOSITORY_DIRECTORY")`, or put
this directory on adopath. Only built-in Stata/Mata required at runtime.

```stata
set type double
insheet using tests/fixtures/data.csv, comma names clear
fasttwfe y x, time(time) cluster(state) gvar(g) inference(native jackknife)
```

Signature: `fasttwfe y [x...] [if] [in] [aw=w], time(t) cluster(c)
gvar(g)|treatment(d) [inference(string) assignments(matrix) rireps(#) seed(#) null(#)]`.
Both timing inputs may be supplied for consistency validation. gvar=0 means never.
Treatment must be cluster-time constant and absorbing. Missing/nonfinite data,
nonpositive weights, rank failures and failed deletions are explicit errors.

FE are cluster and time, using contemporaneous common-slope covariates except
historical FLEX interacted slopes. RCS is supported; individual panel FE/ivar()
is explicitly unsupported in this Stata milestone. Panel data may be analyzed
only as a cluster-FE specification, not relabeled individual FE. Dynamic effects,
event windows, harmonized targets, separate assignment units and alternate
historical options are unsupported.

Use double-precision CSV import: `set type double` before `insheet`.
The test machine's import delimited required unavailable Java; insheet avoids that
development-environment dependency. Tests contain no licensed/proprietary data.
