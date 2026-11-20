version 16.0
clear all
set more off
set type double
capture set processors 1
adopath ++ "."
insheet using "tests/fixtures/data.csv", comma names clear
matrix A=(0,0,0,0,3,3,4,4 \ 4,0,0,0,0,3,3,4 \ 4,4,0,0,0,0,3,3 \ 3,4,4,0,0,0,0,3)
tempname fh M L
file open `fh' using "tests/stata_results.csv", write replace
file write `fh' "cov,metric,index,value" _n
forvalues cov=0/1 {
    local x ""
    if `cov'==1 local x "x"
    fasttwfe y `x', time(time) cluster(state) gvar(g) ///
        inference(native jackknife ri_beta ri_t_native ) assignments(A)
    foreach metric in estimate native_se native_p jackknife_se ri_beta_p ri_t_native_p  {
        file write `fh' "`cov',`metric',0," %24.17g (e(`metric')) _n
    }
    matrix L=e(deletions)
    forvalues j=1/8 {
        file write `fh' "`cov',deleted,`j'," %24.17g (L[1,`j']) _n
    }
    matrix M=e(randomized_statistics)
    forvalues j=1/4 {
        forvalues k=1/2 {
            local idx=(`j'-1)*2+`k'
            file write `fh' "`cov',randomized,`idx'," %24.17g (M[`j',`k']) _n
        }
    }
}
file close `fh'
assert "`e(vce)'"=="cluster"
fasttwfe y x, time(time) cluster(state) gvar(g) inference(jackknife)
assert "`e(vce)'"=="jackknife"
assert missing(e(native_se)) & missing(e(native_p))
assert abs(e(V)[1,1]-e(jackknife_se)^2)<1e-12
capture fasttwfe y, time(time) cluster(state) gvar(g) inference(madeup)
assert _rc==198
preserve
replace d=1 in 1
capture fasttwfe y, time(time) cluster(state) treatment(d)
assert _rc!=0
restore
fasttwfe y x, time(time) cluster(state) treatment(d) inference(point)
assert e(estimate)<.
fasttwfe y, time(time) cluster(state) gvar(g) inference(ri_beta) assignments(A)
assert e(ri_beta_p)>0 & e(ri_beta_p)<=1
matrix badA=A\A[1,1...]
capture fasttwfe y, time(time) cluster(state) gvar(g) inference(ri_beta) assignments(badA)
assert _rc!=0
preserve
drop if state==4
fasttwfe y, time(time) cluster(state) gvar(g) inference(jackknife)
assert e(jackknife_se)>0 & e(jackknife_se)<.
restore

fasttwfe y x [aw=w], time(time) cluster(state) gvar(g) inference(native)
scalar theta=e(estimate)
scalar se=e(native_se)
regress y x d i.state i.time [aw=w], vce(cluster state)
assert abs(theta-_b[d])<1e-9
assert abs(se-_se[d])<1e-9
local initialrng "`c(rngstate)'"
fasttwfe y x, time(time) cluster(state) gvar(g) ///
    inference(ri_beta ri_t_native ) rireps(3) seed(923)
assert "`c(rngstate)'"=="`initialrng'"
assert e(ri_reps)==3
matrix Generated=e(assignments)
scalar saved_ri=e(ri_beta_p)
fasttwfe y x, time(time) cluster(state) gvar(g) ///
    inference(ri_beta ri_t_native ) assignments(Generated)
assert abs(saved_ri-e(ri_beta_p))<1e-12
fasttwfe y x, time(time) cluster(state) gvar(g) ///
    inference(ri_beta ri_t_native ) rireps(3) seed(923)
mata: st_numscalar("draw_difference",max(abs(st_matrix("Generated")-st_matrix("e(assignments)"))))
assert draw_difference==0
mata: st_matrix("Exact",fasttwfe_v1_draws((0\0\3\3),999))
assert rowsof(Exact)==5
mata: Sampled=fasttwfe_v1_draws((J(10,1,0)\J(5,1,3)\J(5,1,4)),4)
mata: st_numscalar("sample_ok",rows(Sampled)==4 & all(rowsum(Sampled:==0):==10) & all(rowsum(Sampled:==3):==5) & all(rowsum(Sampled:==4):==5))
assert sample_ok==1
local initialrng "`c(rngstate)'"
capture fasttwfe y, time(time) cluster(state) gvar(g) inference(ri_beta) assignments(badA) seed(746)
assert _rc!=0
assert "`c(rngstate)'"=="`initialrng'"
display "FOUR_PORT_TEST_EXECUTION_PASS"
