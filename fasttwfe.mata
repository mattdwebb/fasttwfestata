version 16.0
mata:
real scalar fasttwfe_v1_version() return(1)

void fasttwfe_v1_fail(string scalar why)
{
    errprintf("%s\n", why)
    _error(3498)
}

real colvector fasttwfe_v1_codes(real colvector x)
{
    real colvector u, z
    real scalar i
    u=uniqrows(sort(x,1)); z=J(rows(x),1,.)
    for(i=1;i<=rows(u);i++) z[selectindex(x:==u[i])]=J(sum(x:==u[i]),1,i)
    return(z)
}

real matrix fasttwfe_v1_dummies(real colvector x)
{
    real matrix z
    real scalar i
    z=J(rows(x),max(x),0)
    for(i=1;i<=rows(x);i++) z[i,x[i]]=1
    return(z)
}

real colvector fasttwfe_v1_sum(real colvector x, real colvector c)
{
    real colvector z
    real scalar i
    z=J(max(c),1,0)
    for(i=1;i<=rows(x);i++) z[c[i]]=z[c[i]]+x[i]
    return(z)
}

real matrix fasttwfe_v1_absorb(real matrix z, real colvector c, real colvector t,
                       real colvector w)
{
    real matrix old, means
    real colvector g, mass
    real scalar it, j, i
    real rowvector scale
    scale=colmax(abs(z)):+1
    for(it=1;it<=10000;it++) {
        old=z
        for(j=1;j<=2;j++) {
            g=(j==1 ? c : t)
            mass=fasttwfe_v1_sum(w,g); means=J(max(g),cols(z),0)
            for(i=1;i<=rows(z);i++) means[g[i],]=means[g[i],]+w[i]*z[i,]
            means=means:/mass
            z=z-means[g,]
        }
        if(max(abs(z-old):/scale)<1e-12) return(z)
    }
    fasttwfe_v1_fail("NUMERICAL_FAILURE: FE absorption nonconvergence")
    return(z)
}

/* Columns: outcome, time, cluster code, weight, covariates. */
real rowvector fasttwfe_v1_fit(real matrix input, real colvector a, real scalar omit,
                       real scalar native, string scalar family)
{
    real colvector y,w,c,t,g,d,u,e,b,scores,h,direction,ix,cell,cells,effect,l,res,counts,smart
    real matrix x,z,xx,base,v,pi,pp,p
    real scalar n,G,q,rank,theta,se,mass,mult,j,k,off,cc,gg,tt,K,i,root1,root2
    real colvector parent,roots
    p=input
    if(omit>0) {
        p=select(p,p[,3]:!=omit)
    }
    y=p[,1]; w=p[,4]; c=fasttwfe_v1_codes(p[,3]); t=fasttwfe_v1_codes(p[,2]); g=a[p[,3]]
    x=(cols(p)>4 ? p[,5..cols(p)] : J(rows(p),0,.))
    d=(g:>0):&(p[,2]:>=g); u=1:-d
    n=rows(p); G=max(c); q=cols(x)+1
    if(sum(d)==0 | sum(u)==0) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: treatment support")
    mult=1; se=.
    if(family=="twfe") {
        z=fasttwfe_v1_absorb((y,x,d),c,t,w)
        xx=z[,2..cols(z)]
        pi=pinv(xx:*sqrt(w),rank)
        if(rank!=cols(xx)) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: rank deficient")
        b=pi*(sqrt(w):*z[,1]); theta=b[rows(b)]
        if(native) {
            parent=(1::(G+max(t)))
            for(i=1;i<=n;i++) {
                root1=c[i]; root2=G+t[i]
                while(parent[root1]!=root1) root1=parent[root1]
                while(parent[root2]!=root2) root2=parent[root2]
                parent[root2]=root1
            }
            roots=parent
            for(i=1;i<=rows(parent);i++) {
                root1=i
                while(parent[root1]!=root1) root1=parent[root1]
                roots[i]=root1
            }
            K=G+max(t)-rows(uniqrows(sort(roots,1)))+cols(xx)
            if(G<2 | n<=K) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: degrees of freedom")
            e=z[,1]-xx*b
            scores=fasttwfe_v1_sum(sqrt(w):*e:*pi[rows(pi),]',c)
            mult=G/(G-1)*(n-1)/(n-K)
        }
    }
    else {
        K=(family=="flex" ? (G+max(t)-1+sum(d))*q : G+max(t)-1+cols(x))
        if(n*K>12000000) fasttwfe_v1_fail("UNSUPPORTED: reference size cap")
        base=fasttwfe_v1_dummies(c)
        z=fasttwfe_v1_dummies(t)
        if(cols(z)>1) base=(base,z[,2..cols(z)])
        if(family=="did2s" | family=="bjs") {
            z=(base,x); ix=selectindex(u)
            pi=pinv(z[ix,]:*sqrt(w[ix]),rank)
            if(rank!=cols(z)) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: untreated rank")
            b=pi*(sqrt(w[ix]):*y[ix]); e=y-z*b
            mass=sum(w:*d); theta=sum(w:*d:*e)/mass
            if(native) {
                h=quadcolsum(z:*(w:*d))'/mass
                direction=(pi'*h):*sqrt(w[ix])
                if(family=="did2s") res=w:*d:*(e:-theta)/mass
                else {
                    if(any(w:!=w[1])) fasttwfe_v1_fail("UNSUPPORTED: BJS native weights")
                    res=e
                    for(gg=1;gg<=rows(a);gg++) if(a[gg]>0) {
                        for(tt=1;tt<=max(t);tt++) {
                            l=selectindex((g:==a[gg]):&(t:==tt):&d)
                            if(rows(l)>0) {
                                counts=fasttwfe_v1_sum(J(rows(l),1,1),c[l])
                                smart=counts[c[l]]
                                res[l]=e[l]:-sum(smart:*e[l])/sum(smart)
                            }
                        }
                    }
                    res=w:*d:*res/mass
                }
                res[ix]=-direction:*e[ix]
                scores=fasttwfe_v1_sum(res,c)
            }
        }
        else {
            v=(J(n,1,1),x); z=J(n,0,.)
            for(j=1;j<=cols(base);j++) z=(z,base[,j]:*v)
            off=cols(z)
            cell=(c:-1)*max(t)+t
            cells=uniqrows(sort(select(cell,d),1))
            for(j=1;j<=rows(cells);j++) z=(z,(cell:==cells[j]):*v)
            pi=pinv(z:*sqrt(w),rank); b=pi*(sqrt(w):*y)
            h=J(cols(z),1,0); effect=J(n,1,0); mass=sum(w:*d)
            for(j=1;j<=rows(cells);j++) {
                ix=selectindex(cell:==cells[j]); l=(off+(j-1)*q+1::off+j*q)
                h[l]=sum(w[ix])/mass*mean(v[ix,])'
                effect[ix]=v[ix,]*b[l]
            }
            if(max(abs(h'*pi*(z:*sqrt(w))-h'))>1e-8)
                fasttwfe_v1_fail("IDENTIFICATION_FAILURE: FLEX contrast")
            theta=sum(h:*b)
            if(native) {
                if(any(w:!=w[1])) fasttwfe_v1_fail("UNSUPPORTED: FLEX native weights")
                if(G<2 | n<=rank) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: degrees of freedom")
                scores=fasttwfe_v1_sum(sqrt(w):*(y-z*b):*(h'*pi)' + w:*d:*(effect:-theta)/mass,c)
                mult=G/(G-1)*(n-1)/(n-rank)
            }
        }
    }
    if(native) se=sqrt(sum(scores:^2)*mult)
    if(theta>=.) fasttwfe_v1_fail("NUMERICAL_FAILURE: nonfinite estimate")
    return((theta,se))
}

real rowvector fasttwfe_v1_jack(real matrix p, real colvector a, string scalar family)
{
    real colvector dels
    real scalar j,G,se
    real rowvector f
    G=rows(a)
    if(G<2) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: jackknife clusters")
    dels=J(G,1,.)
    for(j=1;j<=G;j++) {
        f=fasttwfe_v1_fit(p,a,j,0,family); dels[j]=f[1]
    }
    se=sqrt((G-1)/G*sum((dels:-mean(dels)):^2))
    if(se<=64*epsilon(1)*max((1\abs(dels))) | se>=.)
        fasttwfe_v1_fail("NUMERICAL_FAILURE: zero jackknife variance")
    return((se,dels'))
}


void fasttwfe_v1_permutations(real rowvector prefix, real colvector left,
                       real rowvector observed, real matrix output)
{
    real colvector values, keep, rest
    real scalar j, at
    if(rows(left)==0) {
        if(any(prefix:!=observed)) output=output\prefix
        return
    }
    values=uniqrows(sort(left,1))
    for(j=1;j<=rows(values);j++) {
        at=selectindex(left:==values[j])[1]
        keep=selectindex((1::rows(left)):!=at)
        rest=(rows(left)==1 ? J(0,1,.) : left[keep])
        fasttwfe_v1_permutations((prefix,values[j]),rest,observed,output)
    }
}

real matrix fasttwfe_v1_draws(real colvector a, real scalar reps)
{
    real colvector values, counts, ord
    real rowvector draw
    real matrix output
    real scalar j,k,space,B,attempt,duplicate
    if(reps<1 | reps!=floor(reps) | reps>100000)
        fasttwfe_v1_fail("INVALID_INPUT: rireps must be integer 1..100000")
    values=uniqrows(sort(a,1)); counts=J(rows(values),1,0)
    for(j=1;j<=rows(values);j++) counts[j]=sum(a:==values[j])
    space=exp(lngamma(rows(a)+1)-sum(lngamma(counts:+1)))
    if(space<1.5) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: constant assignment")
    if(space<=100000.5 & space*rows(a)<=12000000) {
        output=J(0,rows(a),.)
        fasttwfe_v1_permutations(J(1,0,.),a,a',output)
        if(reps>=rows(output)) return(output)
        ord=order(runiform(rows(output),1),1)
        return(output[ord[1..reps],])
    }
    B=min((reps,floor(space+.5)-1))
    output=J(0,rows(a),.)
    for(attempt=1;attempt<=100*B;attempt++) {
        draw=a[order(runiform(rows(a),1),1)]'
        duplicate=all(draw:==a')
        for(k=1;k<=rows(output);k++) if(all(draw:==output[k,])) duplicate=1
        if(!duplicate) output=output\draw
        if(rows(output)==B) return(output)
    }
    fasttwfe_v1_fail("NUMERICAL_FAILURE: unique assignment sampler exhausted")
    return(output)
}
void fasttwfe_v1_run(string scalar vars, string scalar time, string scalar cluster,
             string scalar gvar, string scalar treatment, string scalar weight,
             string scalar touse, string scalar inference, string scalar assignments,
             string scalar out, string scalar delout, string scalar statout,
             string scalar family, real scalar null, real scalar reps,
             string scalar drawout)
{
    real matrix Y,p,A,stats
    real colvector c,t,w,g,d,a,ix,dd
    real scalar G,i,j,native,jk,ri,theta,ns,js,df,tstat,hits,o
    real rowvector f,jack,r
    string rowvector methods
    methods=tokens(inference)
    native=any(methods:=="native") | any(methods:=="ri_t_native")
    jk=any(methods:=="jackknife")
    ri=any(methods:=="ri_beta") | any(methods:=="ri_t_native")
    Y=st_data(.,tokens(vars),touse); t=st_data(.,time,touse)
    c=st_data(.,cluster,touse); w=st_data(.,weight,touse)
    G=max(c); a=J(G,1,0)
    if(gvar!="") g=st_data(.,gvar,touse)
    if(treatment!="") d=st_data(.,treatment,touse)
    if(treatment!="" && any((d:!=0):&(d:!=1))) fasttwfe_v1_fail("INVALID_INPUT: binary treatment")
    if(any(w:<=0)) fasttwfe_v1_fail("INVALID_INPUT: positive weights required")
    for(i=1;i<=G;i++) {
        ix=selectindex(c:==i)
        if(gvar!="") {
            if(min(g[ix])!=max(g[ix]) | min(g[ix])<0) fasttwfe_v1_fail("INVALID_INPUT: gvar")
            a[i]=g[ix[1]]
        }
        else if(sum(d[ix])>0) a[i]=min(select(t[ix],d[ix]))
        dd=(a[i]>0):&(t[ix]:>=a[i])
        if(treatment!="") {
            if(any(dd:!=d[ix])) fasttwfe_v1_fail("INVALID_INPUT: treatment history")
        }
        if(a[i]>0 && (sum(dd)==0 | sum(1:-dd)==0)) fasttwfe_v1_fail("IDENTIFICATION_FAILURE: pre/post")
    }
    p=(Y[,1],t,c,w)
    if(cols(Y)>1) p=(p,Y[,2..cols(Y)])
    f=fasttwfe_v1_fit(p,a,0,native,family); theta=f[1]; ns=f[2]
    r=J(1,12,.); r[1]=theta; r[11]=rows(p); r[12]=G
    if(native) {
        if(ns<=0 | ns>=.) fasttwfe_v1_fail("NUMERICAL_FAILURE: native variance")
        df=(family=="twfe" | family=="flex" ? G-1 : .)
        tstat=(theta-null)/ns
        r[2]=ns; r[3]=tstat; r[4]=df
        r[5]=(df<. ? 2*ttail(df,abs(tstat)) : 2*normal(-abs(tstat)))
    }
    if(jk) {
        jack=fasttwfe_v1_jack(p,a,family); js=jack[1]
        r[6]=js; r[7]=(theta-null)/js; r[8]=2*ttail(G-1,abs(r[7]))
        st_matrix(delout,jack[,2..cols(jack)])
    }
    if(ri) {
        if(null!=0) fasttwfe_v1_fail("UNSUPPORTED: RI nonzero null")
        if(assignments=="") A=fasttwfe_v1_draws(a,reps)
        else A=st_matrix(assignments)
        if(cols(A)!=G | rows(A)<1 | hasmissing(A)) fasttwfe_v1_fail("INVALID_INPUT: assignments")
        stats=J(rows(A),2,.)
        for(i=1;i<=rows(A);i++) {
            if(any(sort(A[i,]',1):!=sort(a,1)) | all(A[i,]':==a))
                fasttwfe_v1_fail("INVALID_INPUT: assignment counts or observed assignment")
            for(j=1;j<i;j++) if(all(A[i,]:==A[j,])) fasttwfe_v1_fail("INVALID_INPUT: duplicate assignment")
            f=fasttwfe_v1_fit(p,A[i,]',0,any(methods:=="ri_t_native"),family)
            stats[i,1]=f[1]
            if(any(methods:=="ri_t_native")) {
                if(f[2]<=0 | f[2]>=.) fasttwfe_v1_fail("NUMERICAL_FAILURE: randomized native variance")
                stats[i,2]=f[1]/f[2]
            }
        }
        for(j=1;j<=2;j++) {
            if(j==1) o=abs(theta)
            if(j==2) o=abs(r[3])
            if(!hasmissing(stats[,j])) {
                hits=sum(abs(stats[,j]):>=o:-64*epsilon(1)*rowmax((J(rows(stats),1,1),abs(stats[,j]),J(rows(stats),1,o))))
                r[8+j]=(1+hits)/(rows(A)+1)
            }
        }
        st_matrix(statout,stats)
        st_matrix(drawout,A)
    }
    st_matrix(out,r)
}
end
