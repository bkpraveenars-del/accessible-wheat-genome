import numpy as np, pandas as pd, nature_fig as nf
from scipy.stats import spearmanr
nf.apply(); W=nf.WONG
d=pd.read_csv("../data/Source_Data_Fig2.csv")
d["delta"]=d.fancm_cM-d.control_cM; d["fc"]=d.fancm_cM/d.control_cM; d["l2"]=np.log2(d.fc)
d["cls"]=np.where(d.fc>1.1,"up",np.where(d.fc<0.9,"down","nc"))
col={"up":W["blue"],"down":W["vermillion"],"nc":nf.GREY}
fig=nf.figure("double",72)
gs=fig.add_gridspec(1,3,left=.065,right=.985,top=.88,bottom=.17,wspace=.42,width_ratios=[1.35,1,1])
a=fig.add_subplot(gs[0]); b=fig.add_subplot(gs[1]); c=fig.add_subplot(gs[2])
# a waterfall of log2 fold change
s=d.sort_values("l2").reset_index(drop=True); x=np.arange(len(s))
a.bar(x,s.l2,color=[col[k] for k in s.cls],width=.75,linewidth=0)
a.axhline(0,color="black",lw=.5)
a.set_xticks(x); a.set_xticklabels(s.interval,rotation=90,fontsize=5)
a.set_ylabel("Map length, $fancm$ / control (log"+r"$_{\mathregular{2}}$"+")")
a.set_xlim(-.7,len(s)-.3)
a.text(.02,.97,"Summed map: %.1f → %.1f cM (+%.1f%%)"%(d.control_cM.sum(),d.fancm_cM.sum(),100*(d.fancm_cM.sum()/d.control_cM.sum()-1)),transform=a.transAxes,va="top",fontsize=6)
for k,lab,yy in [("up","> +10%% (n = %d)"%(s.cls=="up").sum(),.88),("down","< −10%% (n = %d)"%(s.cls=="down").sum(),.81),("nc","within ±10%% (n = %d)"%(s.cls=="nc").sum(),.74)]:
    a.add_patch(__import__("matplotlib").patches.Rectangle((.03,yy-.035),.025,.04,transform=a.transAxes,color=col[k],lw=0))
    a.text(.07,yy-.015,lab,transform=a.transAxes,va="center",fontsize=5.5)
a.set_xlabel("Interval (21 measured regions)")
# b scatter
b.loglog(d.control_cM,d.fancm_cM,"o",ms=3.2,mfc="white",mec=W["blue"],mew=.8)
lim=[2.5,150]; b.plot(lim,lim,color=nf.GREY,lw=.6,ls="--"); b.set_xlim(lim); b.set_ylim(lim)
rho=spearmanr(d.control_cM,d.fancm_cM)[0]
b.text(.04,.96,"Spearman $\\rho$ = %.2f\n$n$ = 21 intervals"%rho,transform=b.transAxes,va="top")
b.text(60,40,"1:1",fontsize=5.5,color="black")
for nm in ["7B_2","1A_1","5B_1","1A_2"]:
    r=d[d.interval==nm].iloc[0]; b.annotate(nm,(r.control_cM,r.fancm_cM),xytext=((-14,3) if nm=="1A_1" else (2,4) if nm=="1A_2" else (3,-6) if nm!="7B_2" else (3,3)),textcoords="offset points",fontsize=5)
from matplotlib.ticker import FixedLocator,NullLocator,FixedFormatter
for ax_ in (b.xaxis,b.yaxis):
    ax_.set_major_locator(FixedLocator([3,10,30,100])); ax_.set_major_formatter(FixedFormatter(["3","10","30","100"])); ax_.set_minor_locator(NullLocator())
b.set_xlabel("Control map length (cM)"); b.set_ylabel("$fancm$ map length (cM)")
# c contribution to net gain
s2=d.sort_values("delta",ascending=False).reset_index(drop=True)
cum=np.cumsum(s2.delta)/s2.delta.sum()*100
c.bar(np.arange(21),s2.delta,color=[col[k] for k in s2.cls],width=.75,lw=0)
c.axhline(0,color="black",lw=.5)
c.set_ylabel("Change in map length (cM)"); c.set_xlabel("Intervals ranked by change")
c.set_xticks([0,5,10,15,20]); c.set_xticklabels(["1","6","11","16","21"])
gp=s2.delta[s2.delta>0].sum(); gn=-s2.delta[s2.delta<0].sum(); net=s2.delta.sum()
c.text(.97,.97,"Net gain %.1f cM\nGross +%.0f / −%.0f cM\nTop 2 (7B_2, 1A_1):\n%.0f%% of net gain"%(net,gp,gn,100*s2.delta[:2].sum()/net),transform=c.transAxes,ha="right",va="top")
nf.label_panels(fig,[a,b,c],dx_mm=-9)
nf.save(fig,"../figures/Fig2")
