# -*- coding: utf-8 -*-
"""
GSE101521 (Pantazatos et al. 2017) - bulk RNA-seq of human DLPFC (BA9).
29 non-psychiatric controls, 9 donors with major depressive disorder who died by
causes other than suicide (MDD), and 21 who died by suicide (MDD-S).

This is the only dataset examined here that separates suicide from depression.
Deposited expression values are converted to counts per million and analysed as
log2(CPM + 1). Groups are compared pairwise with Welch t-tests and Mann-Whitney U
tests, and in a linear model adjusting for age, sex, RIN, brain pH and post-mortem
interval, with Benjamini-Hochberg correction across the 17 genes examined.

Input : data/gse101521/  (see data/README.md)
        data/metadata/gse101521_samples.csv
Output: results/gse101521_expression.csv
        results/gse101521_tests.csv
        results/gse101521_men_only.csv

Reproduces Figure 7 and Supplementary Tables S7-S9.
"""
import itertools, numpy as np, pandas as pd
from scipy import stats

GENES = {
 "CYP46A1":"ENSG00000036530","CYP2D6":"ENSG00000100197","CYP27A1":"ENSG00000135929",
 "CYP7B1":"ENSG00000172817","CH25H":"ENSG00000138135","HMGCR":"ENSG00000112972",
 "APOE":"ENSG00000130203","SREBF2":"ENSG00000198911",
 "RBFOX3":"ENSG00000167281","SNAP25":"ENSG00000132639","SLC17A7":"ENSG00000104888",
 "GAD1":"ENSG00000128683","GFAP":"ENSG00000131095","AQP4":"ENSG00000171885",
 "PLP1":"ENSG00000123560","CSF1R":"ENSG00000182578","P2RY12":"ENSG00000169313",
}
PRIMARY = ["CYP46A1","CYP2D6"]
ORDER   = ["CON","MDD","MDD-S"]
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(_ROOT, "data", "gse101521") + os.sep
_OUT = os.path.join(_ROOT, "results")
os.makedirs(_OUT, exist_ok=True)
_o = lambda f: os.path.join(_OUT, f)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); run = 1.0
    for i in range(n-1, -1, -1):
        run = min(run, p[o[i]]*n/(i+1)); q[o[i]] = run
    return q

def ols_contrast(y, X, idx):
    """Return beta, se, t, p and df for column idx of the design matrix."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    df = len(y) - np.linalg.matrix_rank(X)
    s2 = res @ res / df
    XtXi = np.linalg.pinv(X.T @ X)
    se = np.sqrt(s2 * XtXi[idx, idx])
    t = beta[idx]/se
    return beta[idx], se, t, 2*stats.t.sf(abs(t), df), df

# ---------- data ----------
cnt = pd.read_csv(D+"GSE101521_totalRNA_counts.csv.gz", index_col=0)
cnt.index = [str(i).split('.')[0] for i in cnt.index]
meta = pd.read_csv(os.path.join(_ROOT,"data","metadata","gse101521_samples.csv")).set_index("column")
meta = meta.loc[[c for c in cnt.columns if c in meta.index]]
cnt = cnt[meta.index]
cpm = cnt.divide(cnt.sum(axis=0), axis=1)*1e6
print("Samples: %d   Genow w macierzy: %d" % cnt.shape[::-1])
print(meta.diagnosis.value_counts().reindex(ORDER).to_string(), "\n")

expr = {}
for sym, ens in GENES.items():
    if ens in cpm.index:
        v = cpm.loc[ens]
        if isinstance(v, pd.DataFrame): v = v.iloc[0]
        expr[sym] = np.log2(v.astype(float).values + 1)
E = pd.DataFrame(expr, index=meta.index)

long = E.stack().rename("log2CPM").reset_index()
long.columns = ["column","gene","log2CPM"]
long = long.merge(meta.reset_index()[["column","gsm","diagnosis","age","sex","rin","pmi","ph"]], on="column")
long.to_csv(_o("gse101521_expression.csv"), index=False, encoding="utf-8-sig")

# ---------- covariates ----------
cov = meta.copy()
for c in ["rin","pmi","ph","age"]:
    cov[c] = pd.to_numeric(cov[c], errors="coerce")
    cov[c] = cov[c].fillna(cov[c].median())
cov["male"] = (cov.sex == "Male").astype(float)
for c in ["rin","pmi","ph","age"]:
    cov[c] = (cov[c]-cov[c].mean())/cov[c].std()

# ---------- tests ----------
rows = []
for sym in E.columns:
    y_all = E[sym].values
    for g1, g2 in itertools.combinations(ORDER, 2):
        m1 = (meta.diagnosis==g1).values; m2 = (meta.diagnosis==g2).values
        x, y = y_all[m1], y_all[m2]
        tt = stats.ttest_ind(x, y, equal_var=False)
        nx, ny = len(x), len(y)
        sp = np.sqrt(((nx-1)*x.var(ddof=1)+(ny-1)*y.var(ddof=1))/(nx+ny-2))
        g_h = (y.mean()-x.mean())/sp * (1 - 3/(4*(nx+ny)-9)) if sp>0 else np.nan
        # OLS with covariates, restricted to these two groups
        sub = m1 | m2
        Xc = np.column_stack([np.ones(sub.sum()), m2[sub].astype(float),
                              cov.loc[sub,["age","male","rin","ph","pmi"]].values])
        b, se, t_adj, p_adj, df_adj = ols_contrast(y_all[sub], Xc, 1)
        rows.append(dict(gene=sym, comparison=f"{g2} vs {g1}", n1=nx, n2=ny,
            mean_group1=x.mean(), mean_group2=y.mean(), log2_difference=y.mean()-x.mean(),
            fold_change=2**(y.mean()-x.mean()), hedges_g=g_h,
            p_welch=tt.pvalue, p_MW=stats.mannwhitneyu(x,y).pvalue,
            beta_adj=b, p_adj=p_adj, fc_adj=2**b))
res = pd.DataFrame(rows)

# FDR: osobno w kazdym porownaniu, po genach pierwszorzedowych + reszcie panelu
res["q_welch"] = np.nan; res["q_adj"] = np.nan
for cmp_ in res.comparison.unique():
    m = res.comparison==cmp_
    res.loc[m,"q_welch"] = bh(res.loc[m,"p_welch"].values)
    res.loc[m,"q_adj"]   = bh(res.loc[m,"p_adj"].values)
res.to_csv(_o("gse101521_tests.csv"), index=False, encoding="utf-8-sig")

pd.set_option("display.width", 250, "display.max_columns", 40)
print("="*118); print("PRIMARY TARGETS"); print("="*118)
show = ["gene","comparison","n1","n2","log2_difference","fold_change","hedges_g","p_welch","p_MW","fc_adj","p_adj"]
print(res[res.gen.isin(PRIMARY)][show].to_string(index=False, float_format=lambda v:f"{v:.4g}"))
print("\n"+"="*118); print("RESZTA PANELU - porownanie kluczowe (MDD-S vs CON)"); print("="*118)
print(res[(~res.gen.isin(PRIMARY)) & (res.comparison=="MDD-S vs CON")][show].to_string(index=False, float_format=lambda v:f"{v:.4g}"))

# ---------- sensitivity: men only ----------
print("\n"+"="*118); print("SENSITIVITY: MEN ONLY (as in the RT-qPCR cohort)"); print("="*118)
mm = (meta.sex=="Male").values
r2 = []
for sym in PRIMARY:
    if sym not in E.columns: continue
    y_all = E[sym].values
    for g1,g2 in [("CON","MDD-S"),("CON","MDD"),("MDD","MDD-S")]:
        m1=((meta.diagnosis==g1).values)&mm; m2=((meta.diagnosis==g2).values)&mm
        x,y = y_all[m1], y_all[m2]
        if len(x)<3 or len(y)<3: continue
        tt = stats.ttest_ind(x,y,equal_var=False)
        r2.append(dict(gene=sym, comparison=f"{g2} vs {g1}", n1=len(x), n2=len(y),
                       fold_change=2**(y.mean()-x.mean()), p_welch=tt.pvalue,
                       p_MW=stats.mannwhitneyu(x,y).pvalue))
r2=pd.DataFrame(r2); r2.to_csv(_o("gse101521_men_only.csv"), index=False, encoding="utf-8-sig")
print(r2.to_string(index=False, float_format=lambda v:f"{v:.4g}"))

# ---------- association with RIN and other covariates ----------
print("\n"+"="*118); print("QUALITY CONTROL: association with RIN, pH, PMI and age"); print("="*118)
for sym in PRIMARY:
    if sym not in E.columns: continue
    for c in ["rin","ph","pmi","age"]:
        v = pd.to_numeric(meta[c], errors="coerce")
        ok = v.notna().values
        r,p = stats.pearsonr(E[sym].values[ok], v[ok])
        print(f"  {sym:8s} vs {c:5s}  r = {r:+.3f}   p = {p:.3g}")

print("\nWritten: gse101521_expression.csv, gse101521_tests.csv, gse101521_men_only.csv")
