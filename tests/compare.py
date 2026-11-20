"""Compare actual Stata output against frozen historical Python, not another port."""
import pathlib, csv
import numpy as np
import pandas as pd
root=pathlib.Path(__file__).resolve().parent
got=pd.read_csv(root/"stata_results.csv")
checks=0
maxdiff=0.
for cov in (0,1):
    suffix="_x" if cov else ""
    expected=pd.read_csv(root/"fixtures"/f"expected{suffix}.csv")
    rr=pd.read_csv(root/"fixtures"/f"randomized{suffix}.csv").to_numpy()
    dels=pd.read_csv(root/"fixtures"/f"deletions{suffix}.csv").deleted.to_numpy()
    values=got[got["cov"]==cov]
    for name in expected.columns:
        actual=values[values.metric==name].value.to_numpy()
        diff=np.max(np.abs(actual-expected[name].to_numpy()))
        assert diff<1e-9,(name,cov,diff)
        maxdiff=max(maxdiff,diff);checks+=1
    for name, ref in (("deleted",dels),("randomized",rr.ravel())):
        actual=values[values.metric==name].sort_values("index").value.to_numpy()
        diff=np.max(np.abs(actual-ref))
        assert diff<1e-9,(name,cov,diff)
        maxdiff=max(maxdiff,diff);checks+=1
    theta=expected.estimate[0]
    obs=[theta,theta/expected.native_se[0]]
    for j,name in enumerate(("ri_beta_p","ri_t_native_p")):
        ref=(1+sum(abs(rr[:,j])>=abs(obs[j])-1e-12))/5
        assert abs(values[values.metric==name].value.iloc[0]-ref)<1e-12
        checks+=1
print(f"{checks} frozen-oracle checks PASS; max discrepancy {maxdiff:.3g}")
