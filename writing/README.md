in the .Rprofile, deactivate .libPaths("renv/library/macos/R-4.4/aarch64-apple-darwin20"), and direct the library to the renv.lock file.

# Copyright

© 2026 Lukas Graz, Yannik Schlup
This project is released under the MIT License (see `LICENSE`).

[DOI](https://zenodo.org/records/XXXXXXXX)

10.5281/zenodo.XXXXXXXX

Citation metadata can be found in the CITATION.cff file.

# Setup
- `hypotheses.csv`
  - each Hypothethis == 1 row
  - one variable can hav several rows to do some different plotting/modelling
  - read and preprocessed by `R/read_hyp.R`
- `R/DataPrep.R` Loads data `D` and subsets `is.???` into `R`
- `R/eval_hyp.R`
  - this function is a wrapper for the functrions from: 
    - Model: `R/Models.R` -- here all the modelling is defines as discussted in `notes/Modelling.md`
    - Plot: `R/ggplot_hyp.R`

## Render
In the root of this directory excecute: 

```
quarto render
```

This renders 
- `R/analyse_hypothesis.Qmd`  
  - evaluates hypothethis (plotting, modelling)
  - Saves `cache/AEH.rds` object which can be used to report results in 
- Manuscript -> access it via `firefox ./_manuscript/index.html`
