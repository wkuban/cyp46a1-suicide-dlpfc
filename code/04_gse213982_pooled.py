# -*- coding: utf-8 -*-
"""
GSE213982 (Maitra et al. 2023) - single-nucleus RNA-seq of human DLPFC (BA9).

The deposited matrix contains both the 38 female donors of that study and the 34
male donor labels of GSE144136 processed through the same pipeline. Male labels
were matched to GSE144136 donors by cell-barcode overlap; two labels are separate
libraries from one donor and are merged, and one GSE144136 donor is absent, so the
34 labels resolve to 33 male donors. Together with the 38 female donors this gives
71 donors: 34 controls and 37 suicide decedents.

Cell identifier format: F1.AAACCCACACCTCTGT-1.Mic.Mic1
                        donor . barcode . broad class . cluster

NOTE: the MatrixMarket file is 1-based, so both indices are decremented.

Because an apparent difference in a single gene can arise from a shift in overall
transcript composition, CYP46A1 is additionally expressed per donor as a ratio to
neuronal reference transcripts measured in the same nuclei (see 06_normalisation_control.py).

Input : data/gse213982/  (see data/README.md)
        data/metadata/gse213982_donors.csv
Output: results/gse213982_pseudobulk.csv
        results/gse213982_tests.csv
        a cached accumulator file, so the 361-million-entry matrix is read only once

Reproduces Figure 8A-B and Supplementary Tables S10-S11.
"""
import gzip, re, sys, time, numpy as np, pandas as pd
from scipy import stats

import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(_ROOT, "data", "gse213982") + os.sep
_OUT = os.path.join(_ROOT, "results")
os.makedirs(_OUT, exist_ok=True)
_o = lambda f: os.path.join(_OUT, f)

GENES = ["CYP46A1","CYP2D6","CYP27A1","CYP7B1","CH25H","HMGCR","APOE","SREBF2",
         "POR","FDX1","CYP11A1","CYP17A1","CYP19A1",
         "SLC17A7","GAD1","AQP4","PLP1","CSF1R","P2RY12","RBFOX3","SNAP25"]
MARKERS = {"SLC17A7":"ExN","GAD1":"InN","AQP4":"Ast","PLP1":"Oli","CSF1R":"Mic"}
MIN_NUC = 20

t0=time.time()
genes=[l.strip().strip('"') for l in gzip.open(D+"..._g".replace("..._g",
       "GSE213982_combined_counts_matrix_genes_rows.csv.gz"),'rt')][1:]
cells=[l.strip().strip('"') for l in gzip.open(D+
       "GSE213982_combined_counts_matrix_cells_columns.csv.gz",'rt')][1:]
print("genes %d, cells %d  (%.0fs)"%(len(genes),len(cells),time.time()-t0))

meta=pd.read_csv(os.path.join(_ROOT,"data","metadata","gse213982_donors.csv"))
meta["donor"]=meta.donor.str.replace("M24_2","M24",regex=False)
meta=meta.drop_duplicates("donor")
grp=dict(zip(meta.donor,meta.group)); sex=dict(zip(meta.donor,meta.sex))

def parse(n):
    p=n.split('.')
    d=p[0].replace("M24_2","M24")
    return (d, p[2], p[3]) if len(p)>=4 else (d,"?","?")
print("example:", cells[0], "->", parse(cells[0]))

donor=np.empty(len(cells),object); broad=np.empty(len(cells),object)
for i,c in enumerate(cells):
    d,b,_=parse(c); donor[i]=d; broad[i]=b
uD=sorted(set(donor)); uB=sorted(set(broad))
print("donors %d: %s"%(len(uD),', '.join(uD[:8])+' ...'))
print("classes:", uB)
bad=[d for d in uD if d not in grp]
if bad: sys.exit("MISSING w metadanych: %s"%bad)
di={d:i for i,d in enumerate(uD)}; bi={b:i for i,b in enumerate(uB)}
dIdx=np.array([di[d] for d in donor]); bIdx=np.array([bi[b] for b in broad])
nD,nB=len(uD),len(uB)

gi={}
for g in GENES:
    try: gi[genes.index(g)]=g
    except ValueError: print("  (gene %s not found)" % g)
want=np.zeros(len(genes),bool)
for k in gi: want[k]=True

tot=np.zeros(nD)                      # total counts per donor
totDB=np.zeros((nD,nB))               # total counts per donor x cell class
gD={k:np.zeros(nD) for k in gi}
gDB={k:np.zeros((nD,nB)) for k in gi}
nucDB=np.zeros((nD,nB),int)
for i in range(len(cells)): nucDB[dIdx[i],bIdx[i]]+=1

import os
CACHE=_o("gse213982_accumulators.npz")
if os.path.exists(CACHE):
    _z=np.load(CACHE,allow_pickle=True)
    tot=_z["tot"]; totDB=_z["totDB"]
    gD={k:_z["gD_%d"%k] for k in gi}; gDB={k:_z["gDB_%d"%k] for k in gi}
    print("loaded cached accumulators (%s)"%CACHE)
else:
  print("reading matrix...")
  rd=pd.read_csv(D+"GSE213982_combined_counts_matrix.mtx.gz",sep=r'\s+',skiprows=3,
                 header=None,names=['g','c','v'],dtype={'g':np.int32,'c':np.int32,'v':np.float32},
                 chunksize=12_000_000)
  n=0
  NB=nB
  for ch in rd:
      g=ch['g'].values-1; c=ch['c'].values-1; v=ch['v'].values   # MatrixMarket is 1-based -> convert to 0-based
      d=dIdx[c]; b=bIdx[c]; db=d*NB+b
      tot   += np.bincount(d,  weights=v, minlength=nD)
      totDB += np.bincount(db, weights=v, minlength=nD*NB).reshape(nD,NB)
      m=want[g]
      if m.any():
          gg,dd,ddbb,vv=g[m],d[m],db[m],v[m]
          for k in gi:
              sms=gg==k
              if sms.any():
                  gD[k]  += np.bincount(dd[sms],   weights=vv[sms], minlength=nD)
                  gDB[k] += np.bincount(ddbb[sms], weights=vv[sms], minlength=nD*NB).reshape(nD,NB)
      n+=len(ch)
      if n % 100_000_000 < 20_000_000: print("  %d million (%.0fs)"%(n//1_000_000,time.time()-t0),flush=True)
  print("  read %d million (%.0fs)"%(n//1_000_000,time.time()-t0),flush=True)

  np.savez_compressed(CACHE, tot=tot, totDB=totDB,
      **{"gD_%d"%k:gD[k] for k in gi}, **{"gDB_%d"%k:gDB[k] for k in gi})
  print("accumulators cached -> %s"%CACHE)

rows=[]
for k,name in gi.items():
    for d in range(nD):
        rows.append(dict(gene=name,donor=uD[d],cell_class="ALL",group=grp[uD[d]],
                         sex=sex[uD[d]],n_nuclei=int(nucDB[d].sum()),
                         CPM=gD[k][d]/tot[d]*1e6 if tot[d]>0 else np.nan))
        for b in range(nB):
            if nucDB[d,b]>=MIN_NUC and totDB[d,b]>0:
                rows.append(dict(gene=name,donor=uD[d],cell_class=uB[b],group=grp[uD[d]],
                                 sex=sex[uD[d]],n_nuclei=int(nucDB[d,b]),
                                 CPM=gDB[k][d,b]/totDB[d,b]*1e6))
P=pd.DataFrame(rows); P["log2CPM"]=np.log2(P.CPM+1)
P.to_csv(_o("gse213982_pseudobulk.csv"),index=False,encoding="utf-8-sig")

print("\n=== QC: each marker should be highest in its own cell class ===")
ok=True
for mk,cls in MARKERS.items():
    s=P[(P.gene==mk)&(P.cell_class!="ALL")].groupby("cell_class").log2CPM.mean().sort_values(ascending=False)
    top=s.index[0]; good=top==cls
    ok&=good
    print("  %-8s expected %-4s  highest: %-4s %s"%(mk,cls,top,"OK" if good else "<<< WRONG"))
if not ok: print("  !! cell-class assignment or matrix indexing may be wrong")

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); run=1.0
    for i in range(n-1,-1,-1):
        run=min(run,p[o[i]]*n/(i+1)); q[o[i]]=run
    return q

res=[]
for name in gi.values():
    for cls in ["ALL"]+uB:
        for kohorta,sel in [("women",P.sex=="F"),("men",P.sex=="M"),("pooled",P.sex.notna())]:
            s_=P[(P.gene==name)&(P.cell_class==cls)&sel]
            x=s_.log2CPM[s_.group=="Control"].values; y=s_.log2CPM[s_.group=="Case"].values
            if len(x)<3 or len(y)<3: continue
            tt=stats.ttest_ind(x,y,equal_var=False)
            n1,n2=len(x),len(y)
            sp=np.sqrt(((n1-1)*x.var(ddof=1)+(n2-1)*y.var(ddof=1))/(n1+n2-2))
            g_h=(y.mean()-x.mean())/sp*(1-3/(4*(n1+n2)-9)) if sp>0 else np.nan
            res.append(dict(gene=name,cohort=kohorta,cell_class=cls,n_control=n1,n_case=n2,
                            log2CPM_control=x.mean(),log2CPM_case=y.mean(),
                            fold_change=2**(y.mean()-x.mean()),hedges_g=g_h,
                            t=tt.statistic,df=tt.df,p=tt.pvalue,
                            p_MW=stats.mannwhitneyu(x,y).pvalue))
R=pd.DataFrame(res)
R["q"]=np.nan
for (cls,koh),m in R.groupby(["cell_class","cohort"]).groups.items():
    R.loc[m,"q"]=bh(R.loc[m,"p"].values)
R.sort_values(["gene","cohort","p"]).to_csv(_o("gse213982_tests.csv"),index=False,encoding="utf-8-sig")

pd.set_option("display.width",230,"display.max_columns",30)
show=["gene","cohort","cell_class","n_control","n_case","fold_change","hedges_g","p","p_MW","q"]
print("\n=== PRIMARY TARGETS ===")
print(R[(R.gene.isin(["CYP46A1","CYP2D6","CYP27A1"]))&((R.cell_class=="ALL")|(R.p<0.05))].sort_values(["gene","cohort","p"])[show]
      .to_string(index=False,float_format=lambda v:f"{v:.4g}"))
print("\n=== REMAINDER OF THE PANEL, all nuclei ===")
print(R[(R.cell_class=="ALL")&(~R.gene.isin(["CYP46A1","CYP2D6","CYP27A1"]))].sort_values(["gene","cohort"])[show]
      .to_string(index=False,float_format=lambda v:f"{v:.4g}"))
print("\nnuclei per cell class:")
print(pd.DataFrame(nucDB,index=uD,columns=uB).sum().to_string())
print("\nWritten gse213982_pseudobulk.csv i gse213982_tests.csv  (%.0fs)"%(time.time()-t0))
