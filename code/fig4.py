"""Fig. 4: proposed decision framework (schematic)."""
import nature_fig as nf, matplotlib.patches as mp
nf.apply(); W=nf.WONG
fig=nf.figure("double",112); ax=fig.add_axes([0.01,0.01,0.98,0.96]); ax.axis("off"); ax.set_xlim(0,100); ax.set_ylim(0,62)
def box(x,y,w,h,ec,title,body,tfs=6,bfs=5.5):
    ax.add_patch(mp.Rectangle((x,y),w,h,fc="white",ec=ec,lw=1))
    ax.text(x+w/2,y+h-1.2,title,ha="center",va="top",fontsize=tfs,fontweight="bold")
    ax.text(x+w/2,y+h-4.2,body,ha="center",va="top",fontsize=bfs,linespacing=1.35)
def arr(x1,y1,x2,y2): ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",lw=.7,color="black",mutation_scale=6))
ax.text(0,61,"a",fontsize=8,fontweight="bold",va="top"); ax.text(0,38.5,"b",fontsize=8,fontweight="bold",va="top"); ax.text(0,11.5,"c",fontsize=8,fontweight="bold",va="top")
box(3,45,44,15.5,"black","Record before crossing",
    "Target interval and trait evidence\nParental structure (assembly or breakpoint marker)\nPhase of favourable alleles in the block\nEffective local recombination $r$* in matched material\nFertility or seed-set penalty $w$ of any intervention")
box(53,45,44,15.5,"black","Compute for each candidate cross",
    "Transmissible variance  σ"+nf.sup("2")+"$_{T}$ = $\\mathbf{a}$'$\\mathbf{R}_{\\mathregular{struct}}$$\\mathbf{a}$\nAttainable ceiling ratio  $C$ = $G$*$_{\\mathregular{max}}$ / $G_{\\mathregular{max}}$\nSeeds for target recovery  $N$ = ln β / [$w$ ln(1 − $p$)]\nState the data domain; flag extrapolations", bfs=6)
arr(47,52.75,53,52.75)
ax.plot([75,75],[45,41.5],color="black",lw=.7); ax.plot([11.75,87.75],[41.5,41.5],color="black",lw=.7)
outs=[("Preserve",W["green"],"$C$ ≈ 1; favourable\nalleles in coupling","Transfer the block intact;\nstructural suppression\nprotects the haplotype"),
      ("Break by crossing",W["blue"],"$C$ < 1; $r$* > 0;\n$N$ within budget","Choose a structurally\nmatched parent; size the\npopulation from $N$"),
      ("Change structure",W["orange"],"$C$ < 1; $r$* ≈ 0","Alternative donor, inversion\nreversal, or a modifier only\nif local gain × $w$ > 1"),
      ("Recreate",W["purple"],"Causal change resolved;\nlinked burden high","Edit in the recipient;\nverify expression,\nstructure, phenotype"),
      ("Stop or pre-breed",nf.GREY,"No route meets\n$N$ or $w$ limits","Keep as a pre-breeding\nresource; report the\nnegative result")]
for i,(t,cl,cond,act) in enumerate(outs):
    x=3+i*19.0
    ax.add_patch(mp.Rectangle((x,14.5),17.5,23.5,fc="white",ec=cl,lw=1.2))
    ax.text(x+8.75,36.8,t,ha="center",va="top",fontsize=6,fontweight="bold")
    ax.text(x+8.75,33.2,cond,ha="center",va="top",fontsize=6,linespacing=1.35)
    ax.plot([x+2,x+15.5],[27.0,27.0],color=cl,lw=.6)
    ax.text(x+8.75,25.8,act,ha="center",va="top",fontsize=5.5,linespacing=1.35)
    arr(x+8.75,41.5,x+8.75,38.1)
ax.add_patch(mp.Rectangle((3,1),94,10.5,fc="white",ec="black",lw=.8))
ax.text(50,10.3,"Prospective, falsifiable test",ha="center",va="top",fontsize=6,fontweight="bold")
ax.text(50,7.2,"Compare predicted with observed recovery of the target haplotype, then field gain of recovered lines, against standard genomic cross selection.\n"
        "Reject the framework for a target class if structure-aware predictions do not improve recovery per seed or realised gain in independent crosses.",ha="center",va="top",fontsize=5.5,linespacing=1.4)
for i in range(5): arr(3+i*19+8.75,14.5,3+i*19+8.75,11.6)
nf.save(fig,"../figures/Fig4")
