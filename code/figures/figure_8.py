"""Figure 8: CYP46A1 in the pooled cell-resolved cohort.

GSE213982 together with the GSE144136 male donors it contains (34 controls,
37 suicide decedents). Panel C is the normalisation control: CYP46A1 expressed
per donor as a ratio to SNAP25, SLC17A7 and RBFOX3 measured in the same nuclei.

Inputs : results/gse213982_pseudobulk.csv, results/gse213982_tests.csv
Output : figures/Figure_8 (.tiff/.png/.pdf)
"""
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_RES = os.path.join(_ROOT, "results")
_FIG = os.path.join(_ROOT, "figures")
os.makedirs(_FIG, exist_ok=True)
_r = lambda f: os.path.join(_RES, f)
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family':'Liberation Sans','font.size':9,'axes.linewidth':0.8,
                     'xtick.major.width':0.8,'ytick.major.width':0.8,'pdf.fonttype':42})
CC='#2F6DB5'; SC='#C7782A'; INK='#222222'; MUTED='#6b6b6b'; HI='#8E3B46'

P=pd.read_csv(_r('gse213982_pseudobulk.csv'))
R=pd.read_csv(_r('gse213982_tests.csv'))
ExN=P[P.cell_class=='ExN']
W=ExN.pivot_table(index=['donor','group','sex'],columns='gene',values='log2CPM').reset_index()

fig=plt.figure(figsize=(7.2,6.4))
gs=fig.add_gridspec(2,2,hspace=0.52,wspace=0.42,height_ratios=[1,1])

# ---------- A: CYP46A1 w ExN, wg plci ----------
ax=fig.add_subplot(gs[0,0])
rng=np.random.default_rng(3); pos=0; xt=[]; xl=[]
for sx,lab in [('F','Women'),('M','Men')]:
    for gr,col in [('Control',CC),('Case',SC)]:
        v=W.CYP46A1[(W.sex==sx)&(W.group==gr)].values
        ax.scatter(np.full(len(v),pos)+rng.uniform(-.15,.15,len(v)),v,s=22,
                   facecolor=col,edgecolor='white',linewidth=.5,zorder=3)
        m=v.mean(); h=stats.t.ppf(.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))
        ax.hlines(m,pos-.28,pos+.28,color=INK,lw=1.7,zorder=4)
        ax.errorbar(pos,m,yerr=h,color=INK,capsize=3.5,lw=1,capthick=1,zorder=4)
        xt.append(pos); xl.append('%s\n(n = %d)'%(gr.replace('Case','Suicide'),len(v))); pos+=1
    pos+=0.6
ax.set_xticks(xt); ax.set_xticklabels(xl,fontsize=7.6)
ax.set_ylabel('CYP46A1 in excitatory neurons\nlog$_{\\mathregular{2}}$(CPM + 1)')
ax.spines[['top','right']].set_visible(False)
yl=ax.get_ylim(); ax.set_ylim(yl[0],yl[1]+0.45)
ax.text(0.5,yl[1]+0.28,'Women',ha='center',fontsize=8.5,color=MUTED)
ax.text(3.1,yl[1]+0.28,'Men',ha='center',fontsize=8.5,color=MUTED)
ax.text(-0.26,1.06,'A',transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')
ax.set_title('Pooled: FC 0.80, p = 0.002, q = 0.028',fontsize=8.6,color=INK,pad=3)

# ---------- B: caly panel w ExN, uszeregowany ----------
ax=fig.add_subplot(gs[0,1])
sub=R[(R.cell_class=='ExN')&(R.cohort=='pooled')].copy()
sub['l2']=np.log2(sub.fold_change); sub=sub.sort_values('l2')
mark={'CYP46A1':HI,'SLC17A7':CC,'SNAP25':CC,'RBFOX3':CC}
y=np.arange(len(sub))
cols=[mark.get(g,'#b9b9b9') for g in sub.gene]
ax.axvline(0,color='#9a9a9a',lw=.7,ls=(0,(3,3)),zorder=1)
ax.scatter(sub.l2,y,s=[46 if g in mark else 22 for g in sub.gene],
           c=cols,edgecolor='white',linewidth=.5,zorder=3)
for i,(g,l2) in enumerate(zip(sub.gene,sub.l2)):
    if g in mark:
        ax.annotate(g,(l2,i),xytext=(4 if l2<0 else -4,0),textcoords='offset points',
                    ha='left' if l2<0 else 'right',va='center',fontsize=7.8,
                    color=mark[g],fontweight='bold')
ax.set_yticks([]); ax.set_ylim(-1,len(sub))
ax.set_xlabel('Fold change within excitatory neurons\n(log$_{\\mathregular{2}}$ scale)',fontsize=8.6)
ax.set_xticks([-0.6,-0.3,0,0.3]); ax.set_xticklabels(['0.66','0.81','1','1.23'])
ax.spines[['top','right','left']].set_visible(False)
ax.tick_params(axis='y',length=0)
ax.text(-0.14,1.06,'B',transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')
ax.set_title('All 21 genes assayed, ranked',fontsize=8.6,color=INK,pad=3)

# ---------- C: kontrola normalizacyjna ----------
ax=fig.add_subplot(gs[1,0])
rows=[]
for ref,lab in [(None,'No reference\n(CPM)'),('SNAP25','÷ SNAP25'),('SLC17A7','÷ SLC17A7'),('RBFOX3','÷ RBFOX3')]:
    yv=W.CYP46A1-(W[ref] if ref else 0)
    a=yv[W.group=='Control'].values; b=yv[W.group=='Case'].values
    d=b.mean()-a.mean()
    se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
    df=se**4/((a.var(ddof=1)/len(a))**2/(len(a)-1)+(b.var(ddof=1)/len(b))**2/(len(b)-1))
    tc=stats.t.ppf(.975,df); p=stats.ttest_ind(a,b,equal_var=False).pvalue
    rows.append((lab,d,d-tc*se,d+tc*se,p))
n=len(rows)
ax.axvspan(-0.12,0.12,color='#f0f0f0',zorder=0)
ax.axvline(0,color='#9a9a9a',lw=.8,ls=(0,(3,3)),zorder=1)
for i,(lab,d,lo,hi,p) in enumerate(rows):
    yy=n-1-i; col=HI if p<0.05 else MUTED
    ax.plot([lo,hi],[yy,yy],color=col,lw=1.6,zorder=3,solid_capstyle='butt')
    for xx in (lo,hi): ax.plot([xx,xx],[yy-.12,yy+.12],color=col,lw=1.6,zorder=3)
    ax.plot([d],[yy],marker='o',ms=6.5,color=col,mec='white',mew=.8,zorder=4)
    ax.text(1.03,yy,('p = %.3f'%p) if p>=0.001 else 'p < 0.001',
            transform=ax.get_yaxis_transform(),ha='left',va='center',fontsize=7.8,
            color=col,fontweight='bold' if p<0.05 else 'normal')
ax.set_yticks(range(n)); ax.set_yticklabels([r[0] for r in rows][::-1],fontsize=8.2)
ax.set_ylim(-0.6,n-0.4); ax.set_xlim(-0.92,0.70)
ax.set_xticks([-0.8,-0.4,0,0.4]); ax.set_xticklabels(['0.57','0.76','1','1.32'])
ax.set_xlabel('CYP46A1 ratio, suicide vs control\n(log$_{\\mathregular{2}}$ scale)',fontsize=8.6)
ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
ax.text(-0.42,1.06,'C',transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')

# ---------- D: sklad komorkowy ----------
ax=fig.add_subplot(gs[1,1])
# proportions come from the complete nucleus annotation (08_cell_proportions.py),
# not from the pseudobulk table, which keeps only donor x class with >=20 nuclei
s=pd.read_csv(_r('gse213982_nuclei_counts.csv'))
cls=['ExN','InN','Oli','Ast','OPC','End','Mic','Mix']
w=0.36
for j,(gr,col) in enumerate([('Control',CC),('Case',SC)]):
    m=[s.pct[(s.cell_class==c)&(s.group==gr)].mean() for c in cls]
    e=[s.pct[(s.cell_class==c)&(s.group==gr)].sem() for c in cls]
    ax.bar(np.arange(len(cls))+(j-.5)*w,m,w,yerr=e,color=col,
           error_kw=dict(lw=.8,capsize=2.5),label=gr.replace('Case','Suicide'),zorder=3)
for i,c in enumerate(cls):
    a=s.pct[(s.cell_class==c)&(s.group=='Control')].values; b=s.pct[(s.cell_class==c)&(s.group=='Case')].values
    p=stats.ttest_ind(a,b,equal_var=False).pvalue
    if p<0.05:
        ax.text(i,max(a.mean(),b.mean())+2.0,'*' if p>=0.001 else '**',
                ha='center',fontsize=11,color=INK)
ax.set_xticks(range(len(cls))); ax.set_xticklabels(cls,fontsize=8)
ax.set_ylabel('Share of nuclei per donor (%)',fontsize=8.6)
ax.legend(frameon=False,fontsize=8,loc='upper right')
ax.spines[['top','right']].set_visible(False)
ax.text(-0.26,1.06,'D',transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')

fig.savefig(os.path.join(_FIG,'Figure_8.tiff'), dpi=600, bbox_inches='tight', pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(os.path.join(_FIG,'Figure_8.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(_FIG,'Figure_8.pdf'), bbox_inches='tight')
print('OK')
for r in rows: print('  %-22s FC=%.3f [%.3f, %.3f] p=%.4g'%(r[0].replace('\n',' '),2**r[1],2**r[2],2**r[3],r[4]))
