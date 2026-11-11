"""Fig. 3: analytical model of structurally conditioned within-cross variance. MODEL CALCULATION - no empirical data."""
import numpy as np, nature_fig as nf, matplotlib.patches as mp
from access_model import haldane,R_matrix,var_dh,n_required
nf.apply(); W=nf.WONG
pos=np.arange(20)*2+1.0; blk=(10,30); inb=(pos>blk[0])&(pos<blk[1])
Rc=R_matrix(pos); Rs=R_matrix(pos,blk)
fig=nf.figure("double",118)
gs=fig.add_gridspec(2,4,left=.06,right=.93,top=.93,bottom=.09,hspace=.62,wspace=.55,height_ratios=[.8,1])
a=fig.add_subplot(gs[0,:2]); b1=fig.add_subplot(gs[0,2]); b2=fig.add_subplot(gs[0,3])
c=fig.add_subplot(gs[1,:2]); d=fig.add_subplot(gs[1,2:])
# a: segment schematic, two phase configurations
a.set_xlim(-11,41); a.set_ylim(-.2,3.1); a.axis("off")
cfg={"Coupling block":np.ones(20),"Repulsion block":np.where(inb&(np.arange(20)%2==0),-1,1.0)}
for k,(nm,av) in enumerate(cfg.items()):
    y=2.2-k*1.3
    a.add_patch(mp.FancyBboxPatch((0,y-.09),40,.18,boxstyle="round,pad=0,rounding_size=.09",fc="#D9D9D9",ec="none"))
    a.add_patch(mp.Rectangle((blk[0],y-.13),blk[1]-blk[0],.26,fc="none",ec=W["orange"],lw=1))
    for p,s in zip(pos,av):
        a.plot(p,y+.32,marker="^" if s>0 else "v",ms=3.2,color=W["blue"] if s>0 else W["vermillion"],ls="none")
    ceil_c=np.abs(av).sum(); ceil_s=np.abs(av[~inb]).sum()+abs(av[inb].sum())
    A=var_dh(av,Rs)/var_dh(av,Rc)
    a.text(-10.5,y,nm,va="center",fontsize=6)
    a.text(20,y-.42,"σ"+nf.sup("2")+"$_{T}$/σ"+nf.sup("2")+"$_{\\mathregular{col}}$ = %.3f;  attainable ceiling %.0f%%"%(A,100*ceil_s/ceil_c),ha="center",fontsize=6)
a.add_patch(mp.Rectangle((26,2.93),2,.14,fc="none",ec=W["orange"],lw=1)); a.text(28.6,3.0,"Heterozygous inverted block",va="center",fontsize=5)
a.plot([0,40],[-.1,-.1],color="black",lw=.5); [a.text(x,-.18,str(x),ha="center",va="top",fontsize=5) for x in (0,10,20,30,40)]
a.text(42,-.18,"cM",va="top",fontsize=5)
a.plot([0.5],[3.0],marker="^",color=W["blue"],ms=3,ls="none"); a.text(1.5,3.0,"P1 allele favourable",va="center",fontsize=5)
a.plot([13.5],[3.0],marker="v",color=W["vermillion"],ms=3,ls="none"); a.text(14.5,3.0,"P2 allele favourable",va="center",fontsize=5)
# b: R matrices
for ax,R,t in [(b1,Rc,"Collinear $\\mathbf{R}$"),(b2,Rs,"Structural $\\mathbf{R}_{\\mathregular{struct}}$")]:
    im=ax.imshow(R,cmap="cividis",vmin=0.2,vmax=1,interpolation="nearest")
    ax.set_title(t,fontsize=6); ax.set_xticks([0,19]); ax.set_xticklabels(["1","20"]); ax.set_yticks([0,19]); ax.set_yticklabels(["1","20"])
    ax.set_xlabel("Locus"); ax.set_ylabel("Locus")
cax=fig.add_axes([.945,.62,.008,.28]); cb=fig.colorbar(im,cax=cax); cb.set_ticks([.2,.6,1.0]); cb.set_label("1 − 2$r$*",fontsize=5.5); cb.ax.tick_params(labelsize=5,width=.5); cb.outline.set_linewidth(.5)
# c: random phase configurations
rng=np.random.default_rng(2026); X=[];Y=[];Z=[]
for _ in range(4000):
    av=rng.choice([-1.0,1.0],20)*rng.gamma(2,.5,20)
    X.append(abs(av[inb].sum())/np.abs(av[inb]).sum())
    Y.append(var_dh(av,Rs)/var_dh(av,Rc))
    Z.append((np.abs(av[~inb]).sum()+abs(av[inb].sum()))/np.abs(av).sum())
X,Y,Z=map(np.array,(X,Y,Z))
c.scatter(X,Y,s=1.2,color=W["blue"],alpha=.35,lw=0,label="Within-cross variance ratio")
c.scatter(X,Z,s=1.2,color=W["vermillion"],alpha=.35,lw=0,label="Attainable genetic ceiling (fraction)")
c.axhline(1,color=nf.GREY,lw=.5,ls="--")
c.set_xlabel("Phase coherence of favourable alleles inside block, |Σ$a$|/Σ|$a$|")
c.set_ylabel("Structural / collinear")
c.set_xlim(0,1); c.set_ylim(0,1.6)
c.text(.02,1.52,"Repulsion-rich",fontsize=5.5); c.text(.98,1.52,"Coupling-rich",ha="right",fontsize=5.5)
h=[mp.Patch(color=W["blue"],label="Variance ratio  σ"+nf.sup("2")+"$_{T}$ / σ"+nf.sup("2")+"$_{\\mathregular{col}}$"),mp.Patch(color=W["vermillion"],label="Ceiling ratio  $G$*$_{\\mathregular{max}}$ / $G_{\\mathregular{max}}$")]
c.legend(handles=h,loc="lower right",fontsize=6)
# d: recovery cost
r=np.logspace(-4,-1,200)
for w,ls,lab in [(1,"-","$w$ = 1.00"),(.85,"--","$w$ = 0.85 (fertility penalty)")]:
    d.loglog(r*100,[n_required(x,w) for x in r],color="black",ls=ls,lw=.8,label=lab)
pts=[("1-cM target, collinear",haldane(1),1,W["blue"]),
     ("0.15 cM (YR86 interval, matched parents)",haldane(.15),1,W["green"]),
     ("1 cM with $fancm$ ×2.39 (7B_2 response)",haldane(1)*2.39,.85,W["purple"]),
     ("1 cM with $fancm$ ×0.22 (5B_1 response)",haldane(1)*.22,.85,W["vermillion"])]
for nm,rr,w,cl in pts:
    n=n_required(rr,w); d.plot(rr*100,n,"o",color=cl,ms=3.5); 
d.annotate("",xy=(0.012,2.2e5),xytext=(0.012,4e4),arrowprops=dict(arrowstyle="->",lw=.6))
d.text(0.014,6e4,"$r$* → 0 inside a heterozygous\ninversion: $N$ → ∞",fontsize=5.5)
d.set_xlabel("Effective local recombination $r$* (%)"); d.set_ylabel("DH lines for 95% chance of ≥1 target recombinant")
d.set_xlim(1e-2,10); d.set_ylim(30,3e5)
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator
d.xaxis.set_major_locator(FixedLocator([.01,.1,1,10])); d.xaxis.set_major_formatter(FixedFormatter(["0.01","0.1","1","10"])); d.xaxis.set_minor_locator(NullLocator())
d.yaxis.set_major_locator(FixedLocator([100,1000,1e4,1e5])); d.yaxis.set_major_formatter(FixedFormatter(["100","1,000","10,000","100,000"])); d.yaxis.set_minor_locator(NullLocator())
hh=[__import__("matplotlib").lines.Line2D([],[],color=cl,marker="o",ls="none",ms=3.5,label="%s: %s"%(nm,format(round(n_required(rr,w)),","))) for nm,rr,w,cl in pts]
l1=d.legend(handles=hh,loc="upper right",fontsize=5,handletextpad=.3)
d.add_artist(l1); d.legend(loc="lower left",fontsize=5.5)
nf.label_panels(fig,[a,b1,c,d],dx_mm=-8)
nf.save(fig,"../figures/Fig3")
