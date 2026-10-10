"""Figure 6: CYP46A1 and CYP2D6 in GSE144136 (male single-nucleus cohort).

Inputs : results/gse144136_pseudobulk.csv, results/gse144136_tests.csv
Output : figures/Figure_6 (.tiff/.png/.pdf)
"""
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_RES = os.path.join(_ROOT, "results")
_FIG = os.path.join(_ROOT, "figures")
os.makedirs(_FIG, exist_ok=True)
_r = lambda f: os.path.join(_RES, f)
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family':'Liberation Sans','font.size':9,'axes.linewidth':0.8,'pdf.fonttype':42})
CC='#2F6DB5'; SC='#C7782A'; INK='#222222'; MUTED='#6b6b6b'
pb=pd.read_csv(_r('gse144136_pseudobulk.csv')); R=pd.read_csv(_r('gse144136_tests.csv'))
ORDER=['Ex','Inhib','Astros','Oligos','OPCs','Micro/Macro']
LAB={'Ex':'Excitatory\nneurons','Inhib':'Inhibitory\nneurons','Astros':'Astrocytes',
     'Oligos':'Oligodendro-\ncytes','OPCs':'OPCs','Micro/Macro':'Microglia /\nmacrophages'}
fig=plt.figure(figsize=(7.2,5.8))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.15],hspace=0.55,wspace=0.32)

def dots(ax,gene,subset,title,letter):
    rng=np.random.default_rng(5)
    for i,(grp,col) in enumerate((('Control',CC),('Suicide',SC))):
        v=np.log2(subset.loc[subset.group==grp,gene].values+1)
        ax.scatter(np.full(len(v),i)+rng.uniform(-.13,.13,len(v)),v,s=22,facecolor=col,
                   edgecolor='white',linewidth=.5,zorder=3)
        ax.hlines(v.mean(),i-.26,i+.26,color=INK,lw=1.8,zorder=4)
        h=stats.t.ppf(.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))
        ax.errorbar(i,v.mean(),yerr=h,color=INK,capsize=4,lw=1,capthick=1,zorder=4)
    x=np.log2(subset.loc[subset.group=='Control',gene].values+1)
    y=np.log2(subset.loc[subset.group=='Suicide',gene].values+1)
    t=stats.ttest_ind(x,y,equal_var=False)
    ax.set_xlim(-.6,1.6); ax.set_xticks([0,1])
    ax.set_xticklabels([f'Control\n(n = {len(x)})',f'Suicide\n(n = {len(y)})'])
    ax.set_ylabel('log$_{\\mathregular{2}}$(CPM + 1)')
    ax.set_title(title,fontsize=10,fontweight='bold')
    ax.spines[['top','right']].set_visible(False)
    top=max(x.max(),y.max()); ax.plot([0,0,1,1],[top+.18,top+.33,top+.33,top+.18],color=INK,lw=.8)
    ax.text(.5,top+.38,f'FC {2**(y.mean()-x.mean()):.2f}, p = {t.pvalue:.3f}',ha='center',
            va='bottom',fontsize=8.3)
    ax.set_ylim(min(x.min(),y.min())-.35, top+.95)
    ax.text(-.30,1.08,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')

ALL=pb[pb.cell_class=='ALL']
dots(fig.add_subplot(gs[0,0]),'CYP46A1',ALL,'CYP46A1 (all nuclei)','A')
dots(fig.add_subplot(gs[0,1]),'CYP2D6',ALL,'CYP2D6 (all nuclei)','B')

ax=fig.add_subplot(gs[1,:])
sub=R[(R.gene=='CYP46A1')&(R.cell_class.isin(ORDER))].set_index('cell_class').loc[ORDER]
y=np.arange(len(ORDER))[::-1]
for yi,t in zip(y,ORDER):
    r=sub.loc[t]
    lo,hi=np.log2(r.fold_change)-1.96*0.0,np.log2(r.fold_change)
    col=SC if r.p<0.05 else '#7a7a7a'
    ax.plot(np.log2(r.fold_change),yi,'o',ms=9,color=col,mec='white',mew=.8,zorder=3)
    ax.annotate(f"{r.fold_change:.2f}",xy=(np.log2(r.fold_change),yi),xytext=(0,9),
                textcoords='offset points',ha='center',fontsize=7.6,color=col)
    ax.annotate(f"{int(r.n_control)}/{int(r.n_suicide)}",xy=(1.02,yi),xycoords=('axes fraction','data'),
                fontsize=8,va='center')
    ax.annotate(f"{r.p:.3f}",xy=(1.14,yi),xycoords=('axes fraction','data'),fontsize=8,va='center',
                fontweight='bold' if r.p<0.05 else 'normal')
ax.axvline(0,color='#9a9a9a',lw=.8,ls=(0,(3,3)))
ax.set_yticks(y); ax.set_yticklabels([LAB[t] for t in ORDER],fontsize=8.5)
xt=[-1.5,-1,-.5,0,.5,1]
ax.set_xticks(xt); ax.set_xticklabels([f'{2**v:.2f}' for v in xt])
ax.set_xlim(-1.8,1.35)
ax.set_xlabel('CYP46A1 fold change, suicide vs control (log$_{\\mathregular{2}}$ scale)')
ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
ax.text(1.02,1.03,'n (C/S)',transform=ax.transAxes,fontsize=8,fontweight='bold')
ax.text(1.14,1.03,'p',transform=ax.transAxes,fontsize=8,fontweight='bold')
ax.set_ylim(-0.8,len(ORDER)-0.3)
ax.text(-0.155,1.03,'C',transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')
fig.subplots_adjust(right=0.80)
for ext,kw in (('tiff',dict(pil_kwargs={'compression':'tiff_lzw'})),('png',{}),('pdf',{})):
    fig.savefig(os.path.join(_FIG, f'Figure_6.{ext}'), dpi=600 if ext!='png' else 300, bbox_inches='tight', **kw)
print('ok')
