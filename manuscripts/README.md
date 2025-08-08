## Manuscripts

This folder contains Markdown manuscripts for arXiv/medRxiv/OSF submission:

- `overview_ITC_LLMs.md` — Overview of the ITC advances and LLM-enabled development
- `cNMA_methodology.md` — Component NMA advances
- `MLNMR_methodology.md` — ML‑NMR advances
- `NMI_methodology.md` — NMI advances
- `MAIC_methodology.md` — Bayesian MAIC advances
- `STC_methodology.md` — Bayesian STC advances

### Build PDFs with Pandoc

Requirements: pandoc, LaTeX engine (xelatex recommended), pandoc-crossref, citeproc.

```
cd manuscripts
make
```

PDFs will be generated alongside the sources. Adjust the CSL or template as needed.
