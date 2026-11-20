*! Independent development source 0.1.0. Frozen-oracle reference validation.
program define fasttwfe, eclass sortpreserve
    version 16.0
    syntax varlist(min=1 numeric) [if] [in] [aw], TIME(varname numeric) CLuster(varname) ///
        [GVAR(varname numeric) TREATment(varname numeric) INference(string) ///
        ASSignments(name) NULL(real 0) IVAR(varname) DATAMode(string) ///
        RIReps(integer 0) SEED(integer 1)]
    if "`gvar'"=="" & "`treatment'"=="" {
        display as error "Supply gvar() or treatment()"
        exit 198
    }
    if "`ivar'"!="" {
        display as error "UNSUPPORTED: individual unit FE in this unverified Stata milestone"
        exit 198
    }
    if !inlist("`datamode'","","auto","rcs") {
        display as error "UNSUPPORTED: panel mode requires independently validated unit FE"
        exit 198
    }
    if "`inference'"=="" local inference "native jackknife"
    foreach method of local inference {
        if !inlist("`method'","point","native","jackknife","ri_beta","ri_t_native") {
            display as error "UNSUPPORTED: inference method `method'"
            exit 198
        }
    }
    marksample touse, novarlist
    foreach v of varlist `varlist' `time' `cluster' `gvar' `treatment' {
        quietly count if `touse' & missing(`v')
        if r(N)>0 {
            display as error "INVALID_INPUT: missing observations"
            exit 459
        }
    }
    quietly count if `touse'
    if r(N)==0 error 2000
    tempvar cc w
    quietly egen long `cc'=group(`cluster') if `touse'
    quietly generate double `w'=1 if `touse'
    if "`weight'"!="" quietly replace `w' `exp' if `touse'
    quietly count if `touse' & (missing(`w') | `w'<=0)
    if r(N)>0 error 459
    capture mata: fasttwfe_v1_version()
    if _rc {
        quietly findfile fasttwfe.mata
        quietly run "`r(fn)'"
    }
    tempname result dels stats b V draws
    local savedrng "`c(rngstate)'"
    quietly set seed `seed'
    capture noisily mata: fasttwfe_v1_run("`varlist'","`time'","`cc'","`gvar'","`treatment'","`w'","`touse'","`inference'","`assignments'","`result'","`dels'","`stats'","twfe",`null',`rireps',"`draws'")
    local rc=_rc
    quietly set rngstate `savedrng'
    if `rc' exit `rc'
    matrix colnames `result'=estimate native_se native_stat native_df native_p jackknife_se jackknife_t jackknife_p ri_beta_p ri_t_native_p N N_clust
    matrix `b'=`result'[1,1]
    matrix colnames `b'=effect
    matrix `V'=J(1,1,.)
    matrix colnames `V'=effect
    matrix rownames `V'=effect
    if `result'[1,2]<. matrix `V'[1,1]=`result'[1,2]^2
    else if `result'[1,6]<. matrix `V'[1,1]=`result'[1,6]^2
    if `result'[1,2]<. | `result'[1,6]<. ereturn post `b' `V', esample(`touse')
    else ereturn post `b', esample(`touse')
    if `result'[1,2]<. {
        ereturn local vce "cluster"
        ereturn local vcetype "Cluster"
        if `result'[1,4]<. ereturn scalar df_r=`result'[1,4]
    }
    else if `result'[1,6]<. {
        ereturn local vce "jackknife"
        ereturn local vcetype "Jackknife"
        ereturn scalar df_r=`result'[1,12]-1
    }
    local fields estimate native_se native_stat native_df native_p jackknife_se jackknife_t jackknife_p ri_beta_p ri_t_native_p N N_clust
    local j=0
    foreach f of local fields {
        local ++j
        ereturn scalar `f'=`result'[1,`j']
    }
    if e(jackknife_se)<. ereturn scalar jackknife_df=e(N_clust)-1
    capture confirm matrix `dels'
    if !_rc ereturn matrix deletions=`dels'
    capture confirm matrix `stats'
    if !_rc ereturn matrix randomized_statistics=`stats'
    capture confirm matrix `draws'
    if !_rc {
        ereturn scalar ri_reps=rowsof(`draws')
        ereturn matrix assignments=`draws'
    }
    ereturn local cmd "fasttwfe"
    ereturn local status "PASS"
    ereturn local target_label "TWFE coefficient"
    ereturn local data_mode "rcs"
    display as text "fasttwfe: development reference implementation"
    display as result "estimate = " %12.8f e(estimate)
end
