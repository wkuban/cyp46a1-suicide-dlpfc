"""Figures 2-4: RT-qPCR results.

Figure 2  CYP46A1 and CYP2D6 by group.
Figure 3  Steroidogenic CYPs and the redox partners POR and FDX1.
Figure 4  Summary of all seven estimates with 95% confidence intervals.

Inputs : data/qpcr/prism_raw_values.json, results/qpcr_statistics.json
Outputs: figures/Figure_2, Figure_3, Figure_4 (.tiff/.png/.pdf)
"""
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_RES = os.path.join(_ROOT, "results")
_FIG = os.path.join(_ROOT, "figures")
os.makedirs(_FIG, exist_ok=True)
_r = lambda f: os.path.join(_RES, f)
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy import stats
plt.rcParams.update({'font.family':'Liberation Sans','font.size':9,'axes.linewidth':0.8,'xtick.major.width':0.8,'ytick.major.width':0.8,'pdf.fonttype':42})
d=json.load(open(os.path.join(_ROOT,'data','qpcr','prism_raw_values.json'))); R=json.load(open(_r('qpcr_statistics.json')))
# Short aliases for the keys used below; the JSON uses descriptive names.
_ALIAS = {'fc': 'fold_change', 'fc_lo': 'fold_change_lo', 'fc_hi': 'fold_change_hi',
          'q': 'q_fdr', 'g': 'hedges_g', 'mw_p': 'mannwhitney_p',
          'c_mdd': 'ddct_control_mean', 'c_sd': 'ddct_control_sd',
          's_mdd': 'ddct_suicide_mean', 's_sd': 'ddct_suicide_sd',
          'c_geo': 'geomean_ci_control', 's_geo': 'geomean_ci_suicide'}
for _k in list(R):
    for _short, _long in _ALIAS.items():
        if _long in R[_k]:
            R[_k][_short] = R[_k][_long]
R['ADX'] = R['FDX1']   # the Prism export labels FDX1 as ADX (adrenodoxin)
CC='#2F6DB5'; SC='#C7782A'; INK='#222222'; MUTED='#6b6b6b'
names={'ADX':'FDX1 (adrenodoxin)','POR':'POR'}
def lab2(t): return ('%g'%(2.0**t)) if t>=-3 else '1/%d'%(2**-t)
def save(fig,name):
    fig.savefig(os.path.join(_FIG, name + '.tiff'),dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
    fig.savefig(os.path.join(_FIG, name + '.png'),dpi=300,bbox_inches='tight')
    fig.savefig(os.path.join(_FIG, name + '.pdf'),bbox_inches='tight')
def panel(ax,g,letter):
    c=np.log2(d[g]['raw_ctrl']); s=np.log2(d[g]['raw_sui'])
    rng=np.random.default_rng(7)
    for i,(x,col) in enumerate([(c,CC),(s,SC)]):
        jit=rng.uniform(-0.14,0.14,len(x))
        ax.scatter(np.full(len(x),i)+jit,x,s=24,facecolor=col,edgecolor='white',linewidth=0.6,zorder=3)
        m=x.mean(); h=stats.t.ppf(.975,len(x)-1)*x.std(ddof=1)/np.sqrt(len(x))
        ax.hlines(m,i-0.27,i+0.27,color=INK,lw=1.8,zorder=4)
        ax.errorbar(i,m,yerr=h,color=INK,capsize=4,lw=1,capthick=1,zorder=4)
    allv=np.concatenate([c,s]); lo=np.floor(allv.min())-0.5; hi=np.ceil(allv.max())+1.7
    ax.set_ylim(lo,hi)
    ticks=[t for t in range(int(np.floor(lo)),int(np.ceil(hi))+1)]
    step=2 if len(ticks)>9 else 1
    ticks=[t for t in ticks if t%step==0 and lo<=t<=hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks)); ax.set_yticklabels([lab2(t) for t in ticks])
    ax.axhline(0,color='#9a9a9a',lw=0.6,ls=(0,(3,3)),zorder=1)
    ax.set_xlim(-0.6,1.6); ax.set_xticks([0,1]); ax.set_xticklabels(['Control\n(n = %d)'%len(c),'Suicide\n(n = %d)'%len(s)])
    ax.spines[['top','right']].set_visible(False)
    ax.set_title(names.get(g,g),fontsize=10,fontweight='bold')
    ax.set_ylabel('Relative mRNA expression\n(2$^{\\mathregular{-ΔΔCt}}$, log$_{\\mathregular{2}}$ scale)')
    r=R[g]; ytop=allv.max()+0.75
    ax.plot([0,0,1,1],[ytop-0.2,ytop,ytop,ytop-0.2],color=INK,lw=0.8)
    lab=('p = %.3f'%r['p']) if r['p']>=0.001 else 'p < 0.001'
    ax.text(0.5,ytop+0.12,lab+('*' if r['p']<0.05 else ''),ha='center',va='bottom',fontsize=8.5,color=INK)
    ax.text(-0.34,1.07,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')

# Figure 2
fig,axs=plt.subplots(1,2,figsize=(5.2,3.1))
for ax,g,l in zip(axs,['CYP46A1','CYP2D6'],'AB'): panel(ax,g,l)
fig.tight_layout(w_pad=1.8); save(fig,'Figure_2')
# Figure 3
fig,axs=plt.subplots(2,3,figsize=(7.4,6.2))
for ax,g,l in zip(axs.flat,['CYP11A1','CYP17A1','CYP19A1','POR','ADX'],'ABCDE'): panel(ax,g,l)
axs.flat[5].axis('off')
fig.tight_layout(w_pad=1.6,h_pad=2.2); save(fig,'Figure_3')

# Figure 4 forest plot
order=['CYP46A1','CYP2D6','CYP11A1','CYP17A1','CYP19A1','POR','ADX']
fig,ax=plt.subplots(figsize=(6.6,3.2))
y=np.arange(len(order))[::-1]
for yi,g in zip(y,order):
    r=R[g]; sig=r['p']<0.05
    col=SC if sig else '#7a7a7a'
    ax.plot([np.log2(r['fc_lo']),np.log2(r['fc_hi'])],[yi,yi],color=col,lw=1.6,solid_capstyle='round')
    ax.plot(np.log2(r['fc']),yi,'o',ms=7,color=col,mec='white',mew=0.8,zorder=3)
ax.axvline(0,color='#9a9a9a',lw=0.8,ls=(0,(3,3)))
ax.set_yticks(y); ax.set_yticklabels([names.get(g,g) if g!='ADX' else 'FDX1' for g in order],fontstyle='italic')
xt=[-5,-4,-3,-2,-1,0,1,2]
ax.set_xticks(xt); ax.set_xticklabels([lab2(t) for t in xt]); ax.set_xlim(-5,2.6)
ax.set_xlabel('Fold change, suicide vs control (95% CI, log$_{\\mathregular{2}}$ scale)')
ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
ax.text(1.02,1.04,'n (C/S)',transform=ax.transAxes,fontsize=8,fontweight='bold',ha='left')
ax.text(1.17,1.04,'p',transform=ax.transAxes,fontsize=8,fontweight='bold',ha='left')
ax.text(1.29,1.04,'q (FDR)',transform=ax.transAxes,fontsize=8,fontweight='bold',ha='left')
for yi,g in zip(y,order):
    r=R[g]; yy=(yi+0.5)/len(order)
    tr=ax.get_yaxis_transform()
    ax.annotate('%d/%d'%(r['nc'],r['ns']),xy=(1.02,yi),xycoords=('axes fraction','data'),fontsize=8,va='center')
    ax.annotate(('%.3f'%r['p']),xy=(1.17,yi),xycoords=('axes fraction','data'),fontsize=8,va='center',fontweight='bold' if r['p']<0.05 else 'normal')
    ax.annotate(('%.3f'%r['q']),xy=(1.29,yi),xycoords=('axes fraction','data'),fontsize=8,va='center')
ax.text(-0.15,0.01,'← lower in suicide',transform=ax.get_xaxis_transform(),fontsize=7.5,color=MUTED,ha='right',va='bottom') if False else None
ax.annotate('lower in suicide',xy=(-0.2,-0.95),xycoords=('data','data'),fontsize=7.5,color=MUTED,ha='right')
ax.annotate('higher in suicide',xy=(0.2,-0.95),xycoords=('data','data'),fontsize=7.5,color=MUTED,ha='left')
ax.set_ylim(-1.2,len(order)-0.5)
fig.tight_layout(); save(fig,'Figure_4')

# Figure 5 schematic
fig,ax=plt.subplots(figsize=(8.6,4.7)); ax.set_xlim(0,120); ax.set_ylim(0,62); ax.axis('off')
DOWN=SC; NS='#8a8a8a'; ND='#bdbdbd'
def box(x,y,w,h,text,fc='white',ec='#444444',ls='-',fs=7.6,bold=False,tc=INK):
    b=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.3,rounding_size=1.0',fc=fc,ec=ec,lw=1,ls=ls)
    ax.add_patch(b); ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,fontweight='bold' if bold else 'normal',color=tc)
def arrow(x1,y1,x2,y2,col='#444444'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=9,lw=1.1,color=col,shrinkA=0,shrinkB=0))
def enz(x,y,w,name,detail):
    box(x,y,w,4.2,name,fc='#f6e3cf',ec=DOWN,fs=8.5,bold=True)
    ax.text(x+w/2,y-2.4,detail,ha='center',va='center',fontsize=7,color=DOWN)
ax.text(24,59,'Monoamines & endocannabinoids',ha='center',fontsize=9,fontweight='bold',color=INK)
ax.text(64,59,'Cholesterol elimination',ha='center',fontsize=9,fontweight='bold',color=INK)
ax.text(99,59,'Neurosteroidogenesis',ha='center',fontsize=9,fontweight='bold',color=INK)
for i,(a_,b_) in enumerate([('Tyramine','Dopamine'),('5-Methoxytryptamine','Serotonin'),('Anandamide','Oxidized metabolites')]):
    yy=49-i*7
    box(1,yy,20,4,a_); box(27,yy,20,4,b_); arrow(21.6,yy+2,26.4,yy+2,DOWN)
enz(17,24,14,'CYP2D6','↓ FC 0.25 (95% CI 0.06–0.99), p = 0.049')
ax.plot([24,24],[28.6,34.6],color=DOWN,lw=1,ls=(0,(2,2)))
box(53,49,22,4,'Cholesterol'); box(53,29,22,4,'24S-Hydroxycholesterol'); box(53,18.5,22,4,'Efflux across the BBB',ec='#9a9a9a')
arrow(64,48.4,64,33.6,DOWN); arrow(64,28.4,64,23.1)
enz(66,42,13,'CYP46A1','')
ax.text(74.5,40.2,'↓ FC 0.25 (0.09–0.70)\np = 0.011',ha='center',va='top',fontsize=7,color=DOWN)
ax.text(64,15.2,'CYP46A1 is neuron-specific',ha='center',fontsize=7,color=MUTED,style='italic')
for yy,t in [(49,'Cholesterol'),(39,'Pregnenolone'),(29,'DHEA / androgens'),(19,'Estrogens')]:
    box(90,yy,18,4,t)
for y1,n in [(48.4,'CYP11A1'),(38.4,'CYP17A1'),(28.4,'CYP19A1')]:
    arrow(99,y1,99,y1-5.0)
    ax.text(100.6,y1-2.5,n+'  n.s.',fontsize=7.4,color=MUTED,va='center',fontweight='bold')
ax.text(97.6,45.9,'FDX1 n.s.',fontsize=7,color=MUTED,ha='right',va='center')
box(85.5,7.5,27,4.6,'CYP11B1 · CYP11B2 · CYP21A2',ec=ND,ls='--',fs=6.8)
ax.text(99,4.6,'not detected in DLPFC',ha='center',fontsize=7,color=MUTED)
box(1,8.6,82,4.4,'POR (P450 oxidoreductase), electron donor for CYP2D6, CYP46A1, CYP17A1 and CYP19A1: n.s.',fc='#f4f4f4',ec='#9a9a9a',fs=7.0)
ax.add_patch(FancyBboxPatch((1,0.8),3,2.2,boxstyle='round,pad=0.2',fc='#f6e3cf',ec=DOWN)); ax.text(5.5,1.9,'lower in suicide (nominal p < 0.05)',fontsize=7.2,va='center')
ax.add_patch(FancyBboxPatch((38,0.8),3,2.2,boxstyle='round,pad=0.2',fc='#f4f4f4',ec=NS)); ax.text(42.5,1.9,'no difference (n.s.)',fontsize=7.2,va='center')
ax.add_patch(FancyBboxPatch((62,0.8),3,2.2,boxstyle='round,pad=0.2',fc='white',ec=ND,ls='--')); ax.text(66.5,1.9,'not detected',fontsize=7.2,va='center')
save(fig,'Figure_5')
