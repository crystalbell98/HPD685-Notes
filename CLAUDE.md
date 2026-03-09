# Project: HPD 685 Linear Longitudinal Modeling Notes

## Quarto Rendering
- **Use RStudio's bundled quarto** binary (not system CLI, not `quarto::quarto_render()`):
  ```bash
  /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render
  ```
- Working directory: `/Users/yuan/Documents/HPD 685 Longitudinal/Linear Longitudinal Modeling`

## Deployment
- GitHub Pages: legacy branch mode, serving from `master:/docs`
- Site URL: `https://crystalbell98.github.io/HPD685-Notes/`
- After rendering, commit and push the `docs/` folder

## Project Structure
- `index.qmd` — includes `Longi.md` (longitudinal modeling notes)
- `ols.qmd` — includes `OLS_Regression_Notes.md` (OLS/GLM notes)
- `analysis.qmd` — R code examples (MLM + LSEM measurement invariance)
- `_quarto.yml` — site config (`number-sections: true`, render list)
- `data/` — `antisocial.csv`, `selfconcept.dat`

## Style Conventions
- All content in English
- No numeric prefixes in headings (Quarto adds them via `number-sections: true`)
- No `embed-resources: true` in individual qmd files (website project handles this)
