# Structural variation sets a hidden ceiling on genetic gain in wheat

Code and source data for the Perspective by Praveen Kumar (Department of Genetics and Plant Breeding, College of Agriculture, Agriculture University Jodhpur, India), submitted to *Nature Plants*.

The code reproduces Figs. 1–4 and every numerical value quoted in the main text and Box 1.

## Contents

| Path | Purpose |
|---|---|
| `code/access_model.py` | Box 1 model: structurally conditioned recombination matrix **R**_struct, DH within-cross variance, viable-gamete fraction and recovery cost *N* |
| `code/fig1.py` | Fig. 1, schematic; panel a curves are constrained to published chromosome 3B summary values (Darrier et al. 2017) |
| `code/fig2.py` | Fig. 2, reanalysis of *fancm* interval map lengths |
| `code/fig3.py` | Fig. 3, analytical model figure (no empirical data) |
| `code/fig4.py` | Fig. 4, decision-framework schematic |
| `code/model_values.py` | Prints all model and reanalysis values cited in the text |
| `code/nature_fig.py`, `code/check_nature.py` | Figure style helper and Nature figure-specification checker |
| `code/run_all.py` | Regenerates everything |
| `data/Source_Data_Fig2.csv` | Interval map lengths from Table 1 of Desjardins et al. (2022), *Nat. Commun.* 13, 3644 (CC BY 4.0) |
| `figures/` | Output figures (PDF for submission, SVG with editable text, PNG preview) |

## Reproduce

```bash
pip install -r requirements.txt
cd code
python run_all.py
```

The random configurations in Fig. 3c use a fixed seed (2026), so the output is identical on every run. Fonts fall back to Liberation Sans (metric-identical to Arial) when Arial is not installed.

## Model summary

For DH progeny with additive effects **a**, the within-cross variance is σ²_T = **a**′**R**_struct**a**, with [**R**_struct]_ij = 1 − 2r*_ij. Here r*_ij is computed after setting to zero the map distance inside a heterozygous rearrangement. The attainable ceiling ratio is C = G*_max/G_max. The recovery cost is N = ln β / [w ln(1 − p)]. Full definitions are in Box 1 of the article.

## Data source

Desjardins, S. D. et al. FANCM promotes class I interfering crossovers and suppresses class II non-interfering crossovers in wheat meiosis. *Nat. Commun.* **13**, 3644 (2022). https://doi.org/10.1038/s41467-022-31438-6

## Licence

Code: MIT (see `LICENSE`). Source data in `data/` remain under the CC BY 4.0 licence of the original article.

## Citation

See `CITATION.cff`. The Zenodo concept DOI is [10.5281/zenodo.24134039](https://doi.org/10.5281/zenodo.24134039); the version 1.0.0 DOI is [10.5281/zenodo.24134040](https://doi.org/10.5281/zenodo.24134040).
