"""Structurally conditioned within-cross variance (Box 1). Analytical, no data simulated."""
import numpy as np
def haldane(d_cM): return 0.5*(1-np.exp(-2*np.asarray(d_cM)/100))
def R_matrix(pos, block=None):
    """DH correlation matrix R_ij = 1-2 r*_ij. block=(start,end) cM of a heterozygous inversion:
    crossovers inside it do not yield viable recombinants, so map distance inside the block is removed."""
    pos=np.asarray(pos,float); lo=np.minimum.outer(pos,pos); hi=np.maximum.outer(pos,pos); d=hi-lo
    if block is not None:
        s,e=block; ov=np.clip(np.minimum(hi,e)-np.maximum(lo,s),0,None); d=d-ov
    return 1-2*haldane(d)
def var_dh(a,R): a=np.asarray(a,float); return float(a@R@a)
def w_paracentric(L_cM):
    """Viable-gamete fraction if a single chiasma in an inversion loop unbalances half the gametes (upper-bound loss)."""
    return 1-0.5*(1-np.exp(-2*L_cM/100))
def n_required(r_eff, w=1.0, beta=0.05, target_freq=None):
    """Seeds needed so P(>=1 DH carrying the recombinant coupling haplotype) >= 1-beta. p = r/2 per DH gamete."""
    p = r_eff/2 if target_freq is None else target_freq
    if p<=0: return np.inf
    return np.log(beta)/np.log(1-p)/w
