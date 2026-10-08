"""Figure 7: all available estimates for CYP46A1 and CYP2D6 (forest plot).

Combines this study's RT-qPCR estimates with the re-analyses of GSE144136 and
GSE101521 on a common log2 fold-change scale.

Inputs : results/gse144136_tests.csv, plus the RT-qPCR and GSE101521 estimates
         listed explicitly in this script
Output : figures/Figure_7 (.tiff/.png/.pdf)
"""
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_RES = os.path.join(_ROOT, "results")
_FIG = os.path.join(_ROOT, "figures")
os.makedirs(_FIG, exist_ok=True)
_r = lambda f: os.path.join(_RES, f)
import csv, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family':'Liberation Sans','font.size':9,'axes.linewidth':0.8,
                     'xtick.major.width':0.8,'ytick.major.width':0.8,'pdf.fonttype':42})
QP='#8E3B46'; SN='#2F6DB5'; BK='#C7782A'; INK='#222222'; MUTED='#6b6b6b'

# --- GSE144136: CI z t i df ---
sn={}
for r in csv.DictReader(open(_r('gse144136_tests.csv'))):
    if r['gene'] in ('CYP46A1','CYP2D6'):
        d=np.log2(float(r['fold_change'])); t=float(r['t']); df=float(r['df'])
        se=abs(d/t) if t!=0 else np.nan
        tc=stats.t.ppf(.975,df)
        sn[(r['gene'],r['cell_class'])]=(d,d-tc*se,d+tc*se,int(r['n_control']),int(r['n_suicide']))

L2=lambda fc: np.log2(fc)
ROWS={
'CYP46A1':[
 ('This study — qPCR, bulk DLPFC',            L2(.25), L2(.09), L2(.70), '13/12', QP, 1),
 ('GSE144136 — excitatory neurons',           *sn[('CYP46A1','Ex')][:3],    '17/17', SN, 0),
 ('GSE144136 — inhibitory neurons',           *sn[('CYP46A1','Inhib')][:3], '16/17', SN, 0),
 ('GSE144136 — astrocytes',                   *sn[('CYP46A1','Astros')][:3],'16/12', SN, 0),
 ('GSE144136 — all nuclei (pseudobulk)',      *sn[('CYP46A1','ALL')][:3],   '17/17', SN, 0),
 ('GSE101521 — bulk RNA-seq',                 L2(.923), L2(.762), L2(1.117),'29/21', BK, 0),
],
'CYP2D6':[
 ('This study — qPCR, bulk DLPFC',            L2(.25), L2(.06), L2(.99), '12/12', QP, 1),
 ('GSE144136 — all nuclei (pseudobulk)',      *sn[('CYP2D6','ALL')][:3],   '17/17', SN, 0),
 ('GSE144136 — microglia/macrophages',        *sn[('CYP2D6','Micro/Macro')][:3],'9/13',SN,0),
 ('GSE101521 — bulk RNA-seq',                 L2(1.062), L2(.863), L2(1.308),'29/21', BK, 0),
]}

fig=plt.figure(figsize=(7.1,5.5))
gs=fig.add_gridspec(2,1,height_ratios=[6,4],hspace=0.42)
for k,(gene,letter) in enumerate([('CYP46A1','A'),('CYP2D6','B')]):
    ax=fig.add_subplot(gs[k]); rows=ROWS[gene]; n=len(rows)
    ax.axvspan(-0.3,0.3,color='#f0f0f0',zorder=0)
    ax.axvline(0,color='#9a9a9a',lw=0.8,ls=(0,(3,3)),zorder=1)
    for i,(lab,d,lo,hi,nn,col,diam) in enumerate(rows):
        y=n-1-i
        XL,XR=-3.7,1.0
        clo,chi=max(lo,XL+0.04),min(hi,XR-0.04)
        ax.plot([clo,chi],[y,y],color=col,lw=1.5,zorder=3,solid_capstyle='butt')
        for xx,orig,tip in ((clo,lo,XL),(chi,hi,XR)):
            if abs(orig-xx)>1e-9:   # obcięte – strzałka poza skalę
                ax.annotate('',xy=(tip,y),xytext=(xx,y),zorder=3,
                            arrowprops=dict(arrowstyle='-|>',color=col,lw=1.5,
                                            mutation_scale=9,shrinkA=0,shrinkB=0))
            else:
                ax.plot([xx,xx],[y-.13,y+.13],color=col,lw=1.5,zorder=3)
        if diam:
            ax.plot([d],[y],marker='D',ms=7,color=col,mec='white',mew=.8,zorder=4)
        else:
            ax.plot([d],[y],marker='o',ms=6.5,color=col,mec='white',mew=.8,zorder=4)
        ax.text(1.02,y,nn,transform=ax.get_yaxis_transform(),ha='left',va='center',
                fontsize=8,color=MUTED)
    ax.set_yticks(range(n)); ax.set_yticklabels([r[0] for r in rows][::-1],fontsize=8.5)
    ax.set_ylim(-0.7,n-0.3)
    ax.set_xlim(-3.7,1.0)
    tk=[-3,-2,-1,0,1]
    ax.set_xticks(tk); ax.set_xticklabels(['1/8','1/4','1/2','1','2'])
    ax.set_xlabel('Fold change, suicide vs control (log$_{\\mathregular{2}}$ scale)',fontsize=9)
    ax.set_title(gene,fontsize=10,fontweight='bold',loc='left',pad=6)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.text(-0.52,1.03,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')
    ax.text(1.02,n-0.45,'n (C/S)',transform=ax.get_yaxis_transform(),ha='left',va='center',
            fontsize=8,color=INK,fontweight='bold')

fig.savefig(os.path.join(_FIG,'Figure_7.tiff'), dpi=600, bbox_inches='tight', pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(os.path.join(_FIG,'Figure_7.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(_FIG,'Figure_7.pdf'), bbox_inches='tight')
print('OK')
for g in ROWS:
    for r in ROWS[g]: print('%-8s %-40s FC=%.3f [%.3f, %.3f]'%(g,r[0],2**r[1],2**r[2],2**r[3]))
