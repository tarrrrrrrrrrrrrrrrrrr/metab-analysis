# Statistical Analysis Rules

## Start from the experimental unit

1. Identify factors, levels, blocks/replicates, experimental units, subsamples, seasons/years and sampling times.
2. Treat plot-level replicates as inferential units. Average within-plot leaf, plant or grain subsamples before inferential testing unless a hierarchical model explicitly represents them.
3. Do not infer replication, standard errors or significance from a table of treatment means.
4. Compare the parsed row count with the design equation, for example `2 seasons × 3 cultivars × 3 treatments × 3 field replicates = 54 plots`.

## Model selection

- Factorial randomized complete block: fixed treatment factors plus random block.
- Split plot: main-plot factor tested against main-plot error; subplot factor and interactions tested against subplot error.
- Multi-season experiment: use season as fixed when the stated inference concerns those seasons; use random only when seasons are sampled from a broader target population and the design supports that inference.
- Repeated measurements: represent plot/plant identity and covariance over time. Do not run independent ANOVA at each time without explaining multiplicity.
- Cultivars chosen intentionally are normally fixed effects.

Record the exact mixed model, fixed terms, random terms, denominator degrees-of-freedom method and software versions. Report convergence and singularity warnings.

## Primary inference

For key responses such as yield, 2-AP, Pn, RUE, IPAR and SPAD:

1. Fit the design-correct mixed model.
2. Inspect residuals, leverage, variance patterns and estimability.
3. Report estimated marginal means, standard errors or 95% confidence intervals, exact `n`, effect sizes and adjusted P values.
4. Use Tukey-adjusted compact-letter displays within a clearly defined comparison family. `a/b/c` letters are valid only when calculated from the same fitted model and contrast family.
5. Apply FDR correction when many traits are screened. Keep confirmatory key traits separate from exploratory screening.

## Exploratory analyses

- PCA: standardize only when scale differences justify it; report explained variance and loadings. Do not treat PCA separation as proof of significance.
- Correlation/network/Mantel: disclose pooled factors and treatment confounding; describe association, not causation.
- Regression/path models: report standardization, multicollinearity, model uncertainty and sample size. Do not draw causal arrows from observational covariance alone.
- Response surfaces: show measured-domain limits and observed points. Call a surface descriptive unless validated predictive performance is adequate.
- Random forest/GPR/other ML: use leakage-safe cross-validation; report RMSE, R², folds/repeats, hyperparameters and uncertainty. Negative cross-validated R² means the model is not a strong optimizer.
- GRA/rankings: label as decision-support/exploratory and provide sensitivity checks.

## Prohibited shortcuts

- Pseudoreplication.
- Error bars reconstructed from means.
- Compact-letter display copied from a different model.
- Claiming an optimal range from a weak descriptive surface.
- Calling correlation “activation”, “regulation” or “driving”.
- Altering raw data to match manuscript claims.

## Required comparison table

When reanalysis changes the original paper, create a table with: response, original method, new method, assumption difference, original result, new result, whether the conclusion changes, advantage, limitation and manuscript action.

