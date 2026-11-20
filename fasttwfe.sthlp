{smcl}
{title:fasttwfe: difference-in-differences estimation}
{p}PARTIAL: frozen historical comparisons pass; full official certification pending.{p_end}
{p}{cmd:fasttwfe y [x...] [aw=w], time(t) cluster(c) gvar(g)}
{cmd:inference(native jackknife)}{p_end}
{p}Alternatively supply treatment(d). Optional assignments(matrix), null(#).
Matrix rows are unique nonobserved cohort-label permutations, columns sorted by cluster.
Inference: point native jackknife ri_beta ri_t_native  .
Use rireps(#) seed(#) to generate unique nonobserved cohort-preserving assignments.
Small supports are enumerated; caller RNG state is restored.
Individual panel FE and nonzero RI null are unsupported.
Missing data and unidentified models fail. See README and VALIDATION.{p_end}
