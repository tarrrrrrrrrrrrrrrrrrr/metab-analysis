# Publication Figure Rules

Use Python with matplotlib/seaborn and follow the `nature-figure` skill. Define the conclusion, evidence logic, export needs and reviewer risks before plotting.

## Style contract

- Times New Roman for all figure text where available.
- Use the supplied palette-reference image. Extract its colors once and reuse exact hexadecimal values. If colors are insufficient, extend within the same hue/lightness family and document the additions.
- Use black edges on scatter points, bars, boxplots and filled marks when it improves separation.
- Keep a coherent multi-hue palette; avoid a one-color interface.
- Use real error bars from replicate-level data.
- Add model-supported `a/b/c` compact-letter display labels in empty space above error bars.
- Horizontal tick labels are preferred when they fit. Wrap or abbreviate before rotating.

## Geometry and typography

- Choose standard single-column, 1.5-column or double-column dimensions before plotting.
- Keep panel labels bold and consistently positioned.
- Use stable axes, margins and panel ratios. Do not leave a large empty center between panels.
- Keep legends outside data-dense areas. Never cover bars, points or confidence intervals.
- Avoid nested decorative cards, excessive titles and explanatory text inside figures.

## Chart-specific checks

- Box/violin plots: use wider boxes when readable, show median/IQR clearly, and outline distributions.
- Bars: black edge, mean ± SE/CI, compact letters, and no redundant numerical labels unless requested.
- Heatmaps: cell borders and overlays must match the exact cell grid; top and side labels must have equal spacing.
- Networks/path diagrams: coefficients belong on their corresponding edge. Use a short accurate leader only when the label cannot sit on the edge. Never leave disconnected leaders.
- Response surfaces/3D: show all observed points with black outlines, keep points inside axes limits, and use different markers only for real grouping variables.
- Correlation panels: disclose pooled treatment structure and avoid causal titles.

## Export contract

For every figure:

1. Save each panel separately.
2. Save the combined figure.
3. Export vector PDF and 600 dpi PNG; add SVG when editors need editable vectors.
4. Preserve transparent backgrounds only when intended.
5. Record figure size, DPI, font, palette, error-bar source, significance method and source table in `figure_metadata.csv`.

## QA gate

Check file readability, blank canvases, DPI, dimensions, clipping, missing glyphs, overlapping labels, inconsistent palettes, unsupported significance labels, missing panels, duplicated filenames and manuscript-caption consistency. Inspect rendered output visually before replacing an existing figure.

