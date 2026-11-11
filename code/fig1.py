"""Fig. 1: schematic. Panel a curves are constrained to published summary values (Darrier et al. 2017), not data."""
import numpy as np, nature_fig as nf, matplotlib.patches as mp
from scipy.optimize import brentq
nf.apply(); W=nf.WONG; G=nf.GREY
x=np.linspace(0,1,2001); u=np.abs(x-.5)*2; dx=x[1]-x[0]
dist=u>0.81   # distal 9.5% at each end = 19% of length
def co(k): f=np.exp(k*u); return f/f.sum()/dx
k=brentq(lambda k:(co(k)[dist].sum()*dx)-0.82,1,60)
def gn(c): g=1+c*u**4; return g/g.sum()/dx
c=brentq(lambda c:(gn(c)[~dist].sum()*dx)-0.65,0.01,100)
fig=nf.figure("double",128)
gs=fig.add_gridspec(2,2,left=.07,right=.985,top=.95,bottom=.05,hspace=.38,wspace=.18,height_ratios=[1,1.05])
a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]); cc=fig.add_subplot(gs[1,:])
# a
a.fill_between(x,0,co(k),color=W["blue"],alpha=.55,lw=0,label="Crossovers (82% in distal 19%)")
a.plot(x,gn(c),color=W["orange"],lw=1,label="Genes (≈65% proximal)")
a.axvspan(.095,.905,ymin=0,ymax=1,color="#EFEFEF",lw=0,zorder=0)
a.set_xlim(0,1); a.set_ylim(0,max(co(k))*1.05)
a.set_xticks([0,.095,.5,.905,1]); a.set_xticklabels(["Tel","","Cen","","Tel"])
a.set_xlabel("Relative position along chromosome"); a.set_ylabel("Relative density")
a.text(.5,max(co(k))*.62,"Proximal 81% of length\n≈18% of crossovers,\n≈60–70% of genes",ha="center",fontsize=5.5)
a.legend(loc="upper center",fontsize=5.5,bbox_to_anchor=(.5,1.0))
# b YR86
b.axis("off"); b.set_xlim(0,10); b.set_ylim(0,10)
def chrom(y,segs,lab):
    b.add_patch(mp.FancyBboxPatch((.6,y-.18),7.9,.36,boxstyle="round,pad=0,rounding_size=.18",fc="#D9D9D9",ec="none"))
    b.add_patch(mp.Rectangle((3.4,y-.24),3.0,.48,fc="white",ec=W["orange"],lw=.9))
    for xx in (3.9,4.9,5.9):
        b.annotate("",xy=(xx+(.35 if segs>0 else -.35),y),xytext=(xx-(.35 if segs>0 else -.35),y),arrowprops=dict(arrowstyle="-|>",lw=.7,color="black",mutation_scale=5))
    b.text(8.7,y,lab,va="center",fontsize=5)
b.text(0,9.6,"Structural heterozygote",fontsize=6,fontweight="bold")
chrom(8.7,1,"InvCS"); chrom(7.9,-1,"InvAK58")
b.text(.6,7.05,"Zhongmai 895 × Yangmai 16, Emai 580 or Avocet S",fontsize=5.5)
b.text(.6,6.4,"Crossovers in the 4.1-Mb inverted segment are suppressed or\nyield unbalanced gametes → no usable recombinants",fontsize=5.5,va="top")
b.text(0,4.6,"Structurally matched cross",fontsize=6,fontweight="bold")
chrom(3.7,1,"InvCS"); chrom(2.9,1,"InvCS")
b.plot([5.15,5.65],[3.7,2.9],color=W["vermillion"],lw=1.2); b.plot([5.15,5.65],[2.9,3.7],color=W["vermillion"],lw=1.2)
b.text(.6,2.05,"Zhongmai 895 × Jimai 22",fontsize=5.5)
b.text(.6,1.4,"Normal recombination; $YR86$ mapped to 0.15 cM\n(710.27–712.56 Mb, Chinese Spring RefSeq)",fontsize=5.5,va="top")
# c polyploid
cc.axis("off"); cc.set_xlim(0,20); cc.set_ylim(0,6.2)
cc.text(0,5.95,"Homoeologous buffering",fontsize=6,fontweight="bold")
for i,(g,lab) in enumerate([("A","A"),("B","B"),("D","D")]):
    y=5.1-i*.75
    cc.add_patch(mp.FancyBboxPatch((.4,y-.15),5.6,.3,boxstyle="round,pad=0,rounding_size=.15",fc="#D9D9D9",ec="none"))
    cc.text(.1,y,lab,va="center",ha="center",fontsize=6)
    cc.add_patch(mp.Rectangle((3.0,y-.15),.25,.3,fc=W["blue"],ec="none"))
cc.add_patch(mp.Rectangle((4.3,5.1-1.5-.16),.9,.32,fc="white",ec="black",lw=.5,ls=(0,(2,1))))
cc.text(4.75,5.1-1.5-.45,"deletion",ha="center",va="top",fontsize=5)
cc.text(.4,2.75,"Dosage from A and B tolerates structural change\n(wide design space for editing and introgression)",fontsize=5.5,va="top")
cc.text(.4,1.65,"A recessive allele is masked unless all homoeologs\nare stacked: $P$(access) = $P_{A}$ × $P_{B}$ × $P_{D}$",fontsize=6,va="top")
cc.plot([3.125,3.125],[5.35,3.4],color=W["blue"],lw=.6,ls=":")
cc.text(3.3,5.45,"homoeologous gene copies",fontsize=5)
# gene-pool gradient
cc.text(7.2,5.95,"Gene-pool accessibility for the D and AB subgenomes",fontsize=6,fontweight="bold")
pools=[("$Ae. tauschii$","D-genome donor;\nhomologous with D","Synthetic hexaploids,\ndirect crosses"),
       ("Tetraploid AB relatives","Homologous with\nA and B only","Bridge crosses"),
       ("Wild relatives\n(e.g. $Ae. mutica$)","Homoeologous only;\nsuppressed by\n$Ph1$/$Ph2$","$zip4$-$5B$, $msh7$-$3D$\nor $Ph1$-suppressor\ndonors"),
       ("Engineered structure","Inversion reversal,\ncentromere editing,\nlarge-fragment edits","Verify fertility\nand transmission")]
for j,(t,d1,d2) in enumerate(pools):
    x0=7.2+j*3.2
    cc.add_patch(mp.Rectangle((x0,1.1),2.95,4.25,fc="white",ec=[W["green"],W["skyblue"],W["orange"],W["purple"]][j],lw=1))
    cc.text(x0+1.475,4.95,t,ha="center",va="top",fontsize=5.5,fontweight="bold")
    cc.text(x0+1.475,4.0,d1,ha="center",va="top",fontsize=5.5)
    cc.text(x0+1.475,2.5,d2,ha="center",va="top",fontsize=5.5)
cc.annotate("",xy=(19.95,.6),xytext=(7.2,.6),arrowprops=dict(arrowstyle="-|>",lw=.8,color="black",mutation_scale=7))
cc.text(13.6,.2,"Native recombination with wheat chromosomes decreases; dependence on meiotic or structural intervention increases",ha="center",fontsize=5.5)
nf.label_panels(fig,[a,b,cc],dx_mm=-9)
nf.save(fig,"../figures/Fig1")
print("k",k,"c",c)
