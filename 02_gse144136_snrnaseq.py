import gzip, csv, re, json, time, collections
import numpy as np, pandas as pd
from scipy import stats

import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(_ROOT, "data", "gse144136") + os.sep
_OUT = os.path.join(_ROOT, "results")
os.makedirs(_OUT, exist_ok=True)
_o = lambda f: os.path.join(_OUT, f)

GENES = ["CYP46A1","CYP2D6","CYP27A1","CYP11A1","CYP17A1","CYP19A1","POR","FDX1",
         "P2RY12","CX3CR1","GPR34","SOX9","PMP2",
         "SNAP25","RBFOX3","SLC17A7","GAD1","AQP4","PLP1","CSF1R"]
t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]",*a,flush=True)

def names(p):
    out=[]
    with gzip.open(D+p,'rt') as f:
        r=csv.reader(f); next(r)
        for row in r: out.append(row[-1])
    return out

cells=names("GSE144136_CellNames.csv.gz")
genes=names("GSE144136_GeneNames.csv.gz")
nC,nG=len(cells),len(genes)
log(f"{nC:,} nuclei, {nG:,} genes")

def parse(n):
    p=n.split('_'); cond=p[-3]; batch=p[-2]
    full='_'.join(p[:-3]); cluster,_,donor=full.rpartition('.')
    broad=re.match(r'[A-Za-z/]+',cluster).group(0)
    return broad,cluster,int(donor),cond,batch
meta=[parse(c) for c in cells]
broad=np.array([m[0] for m in meta]); donor=np.array([m[2] for m in meta])
cond=np.array([m[3] for m in meta]); batch=np.array([m[4] for m in meta])

idx={g:i for i,g in enumerate(genes)}
want={idx[g]:g for g in GENES if g in idx}
log("genes found:", sorted(want.values()))
log("MISSING:", [g for g in GENES if g not in idx])
gids=np.array(sorted(want), dtype=np.int32)

total=np.zeros(nC); picks={gi:np.zeros(nC) for gi in want}
n=0
reader=pd.read_csv(D+"GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz", sep=' ',
                   skiprows=2, header=None, names=['g','c','v'],
                   dtype={'g':np.int32,'c':np.int32,'v':np.float32}, chunksize=10_000_000)
for ch in reader:
    g=ch['g'].values-1; c=ch['c'].values-1; v=ch['v'].values
    total += np.bincount(c, weights=v, minlength=nC)
    m=np.isin(g,gids)
    if m.any():
        gm=g[m]; cm=c[m]; vm=v[m]
        for gi in want:
            s=gm==gi
            if s.any(): picks[gi]+=np.bincount(cm[s],weights=vm[s],minlength=nC)
    n+=len(ch); log(f"  {n:,} entries")
log("matrix read")

rows=[]
for key_name, keys in (("by_class", list(zip(donor,broad))), ("all", list(zip(donor,['ALL']*nC)))):
    gr=collections.defaultdict(list)
    for i,k in enumerate(keys): gr[k].append(i)
    for (d,b),ii in gr.items():
        ii=np.array(ii); tot=total[ii].sum()
        if tot<=0: continue
        r={"donor":int(d),"cell_class":b,"group":cond[ii[0]],"batch":batch[ii[0]],
           "n_nuclei":len(ii),"total_counts":float(tot)}
        for gi,sym in want.items(): r[sym]=float(picks[gi][ii].sum()/tot*1e6)
        rows.append(r)
pb=pd.DataFrame(rows)
pb.to_csv(_o("gse144136_pseudobulk.csv"),index=False)
log("pseudobulk written:", pb.shape)

res=[]
for sym in sorted(want.values()):
    for b in sorted(pb.cell_class.unique()):
        sub=pb[(pb.cell_class==b)&(pb.n_nuclei>=20)]
        x=np.log2(sub.loc[sub.group=="Control",sym].values+1)
        y=np.log2(sub.loc[sub.group=="Suicide",sym].values+1)
        if len(x)<5 or len(y)<5: continue
        t=stats.ttest_ind(x,y,equal_var=False)
        res.append({"gene":sym,"cell_class":b,"n_control":len(x),"n_suicide":len(y),
                    "log2CPM_control":x.mean(),"log2CPM_suicide":y.mean(),
                    "fold_change":2**(y.mean()-x.mean()),
                    "t":t.statistic,"df":t.df,"p":t.pvalue,
                    "p_MW":stats.mannwhitneyu(x,y).pvalue})
R=pd.DataFrame(res).sort_values(["gene","p"])
R.to_csv(_o("gse144136_tests.csv"),index=False)
log("tests written:", R.shape)
print(R[R.gene.isin(["CYP46A1","CYP2D6","CYP27A1"])].to_string(index=False))
