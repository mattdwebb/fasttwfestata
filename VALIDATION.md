# Partial validation

Fresh execution under StataNow 19.5, single hidden batch session per package.
One synthetic 8-cluster 5-period unequal-cell fixture, with/without time-varying X.
18 frozen-oracle checks compare point, native SE/p, JK SE, every deletion,
4 supplied assignment statistics and all RI p-values to historical Python.
Source HEAD: 29c4d91805d6024a561ef44be20ca72836464922. Absolute tolerance 1e-9.
Additional Stata assertions: point/RI-only posting, invalid inference/history,
duplicate assignments, cohort-deletion failures, positive-weight points. TWFE
checks weighted point/CR1 against built-in regress; did2s checks weighted first
stage against built-in regress/predict. BJS/FLEX reject native nonuniform weights.