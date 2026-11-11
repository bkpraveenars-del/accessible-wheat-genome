"""Reproduce every numerical value quoted in the main text and Box 1."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from access_model import haldane, R_matrix, var_dh, n_required

d = pd.read_csv("../data/Source_Data_Fig2.csv")
d["delta"] = d.fancm_cM - d.control_cM
top2 = d.nlargest(2, "delta")
print("fancm (Fig. 2, Desjardins et al. 2022 Table 1)")
print(f"  summed map: {d.control_cM.sum():.2f} -> {d.fancm_cM.sum():.2f} cM (+{100*(d.fancm_cM.sum()/d.control_cM.sum()-1):.1f}%)")
print(f"  net gain {d.delta.sum():.2f} cM; top two {top2.interval.tolist()} = {100*top2.delta.sum()/d.delta.sum():.1f}% of net")
print(f"  Spearman rho control vs fancm = {spearmanr(d.control_cM, d.fancm_cM)[0]:.2f}")
print(f"  fold change 7B_2 = {70.99/29.68:.2f}, 5B_1 = {4.51/20.53:.2f}, 1A_2 = {11.44/3.33:.2f}")

pos = np.arange(20) * 2 + 1.0; blk = (10, 30); inb = (pos > blk[0]) & (pos < blk[1])
Rc, Rs = R_matrix(pos), R_matrix(pos, blk)
print("\nModel (Fig. 3, Box 1)")
for name, a in [("coupling", np.ones(20)), ("repulsion", np.where(inb & (np.arange(20) % 2 == 0), -1, 1.0))]:
    vc, vs = var_dh(a, Rc), var_dh(a, Rs)
    rc, rs = (1 - Rc) / 2, (1 - Rs) / 2; off = ~np.eye(20, dtype=bool)
    dfor = 2 * np.sum((np.outer(a, a) * (rc - rs))[off])
    ceil = (np.abs(a[~inb]).sum() + abs(a[inb].sum())) / np.abs(a).sum()
    print(f"  {name}: variance ratio {vs/vc:.3f}; ceiling {100*ceil:.0f}%; Box 1 formula check {dfor:.6f} vs {vs-vc:.6f}")
rng = np.random.default_rng(2026); X, Y, Z = [], [], []
for _ in range(4000):
    a = rng.choice([-1.0, 1.0], 20) * rng.gamma(2, .5, 20)
    X.append(abs(a[inb].sum()) / np.abs(a[inb]).sum()); Y.append(var_dh(a, Rs) / var_dh(a, Rc))
    Z.append((np.abs(a[~inb]).sum() + abs(a[inb].sum())) / np.abs(a).sum())
X, Y, Z = map(np.array, (X, Y, Z))
print(f"  4,000 configurations: rho(coherence, ceiling) = {spearmanr(X, Z)[0]:.2f}; rho(coherence, variance) = {spearmanr(X, Y)[0]:.2f}")
print(f"  variance not reduced by >10% in {100*np.mean(Y >= 0.9):.0f}%; of those, ceiling <80% in {100*np.mean(Z[Y >= 0.9] < 0.8):.0f}%")
print("\nRecovery cost N (95% chance of >=1 target recombinant, DH)")
print(f"  1 cM collinear: {n_required(haldane(1)):.0f}")
print(f"  0.15 cM: {n_required(haldane(.15)):.0f}")
print(f"  1 cM x 2.39 (7B_2), w = 0.85: {n_required(haldane(1)*70.99/29.68, .85):.0f}")
print(f"  1 cM x 0.22 (5B_1), w = 0.85: {n_required(haldane(1)*4.51/20.53, .85):.0f}")
