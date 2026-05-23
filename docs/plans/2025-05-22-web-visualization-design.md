# Cerebral Hemodynamic Simulator — Web Visualization Design

## Framing

Bridge between academic research and modern web engineering. The science is front and center: proper anatomical labels, physiological units, and enough medical context for a technical viewer to understand the thesis contribution.

## Data Pipeline

Pre-compute a grid of simulation scenarios:
- 2 conditions (NSR, AF) x 6 Circle of Willis variants x 3 heart rates (60, 75, 90) = 36 runs
- From each run (~500K timesteps), extract 10 steady-state beats (beats 100-109, after transients settle) — ~800 timesteps per scenario
- Output columns: time, P_a, q_ml, q_al, q_pl, q_mr, q_ar, q_pr
- Save as single JSON keyed by `{condition}_{cow}_{hr}` to `web/public/data/scenarios.json`
- Orchestrated by `scripts/precompute.py`
- Estimated payload: ~2 MB uncompressed, ~200 KB gzipped

## UI Layout

Single-page app, three panels:

```
+-----------------------------------------------------+
|  Cerebral Hemodynamic Simulator                      |
|  Controls: Condition, CoW variant, HR, Compare mode  |
+-------------------------+---------------------------+
|  Aortic Pressure        |  Circle of Willis         |
|  (D3 line chart)        |  (SVG anatomical diagram) |
|  ~10 beats, mmHg        |  flow-colored vessels     |
+-------------------------+---------------------------+
|  Cerebral Blood Flows (6 territories, 2x3 grid)     |
|  MCA-L  ACA-L  PCA-L                                |
|  MCA-R  ACA-R  PCA-R                                |
+-----------------------------------------------------+
```

### Controls
- Two scenario selectors (left/right) for compare mode
- Dropdowns: Condition (NSR/AF), CoW variant (clinical names), Heart Rate (60/75/90)
- Compare toggle: overlays two scenarios with different line styles
- All switching is instant (pre-computed data)

### Aortic Pressure Panel (top-left)
- D3 line chart, ~10 beats
- X-axis: seconds, Y-axis: mmHg
- Compare mode: two traces overlaid (solid vs dashed)

### Circle of Willis Panel (top-right)
- SVG schematic: ICA L/R, BA, MCA L/R, ACA1/2 L/R, PCA1/2 L/R, PCoA L/R, ACoA
- Vessel stroke width scales with mean flow magnitude
- Color encodes magnitude (diverging scale)
- Missing vessels (per CoW variant) shown as dashed grey
- Hover tooltip: vessel name + mean flow in mL/s

### CoW Variant Mapping
| Code | Clinical Name | Missing Vessel(s) |
|------|--------------|-------------------|
| 0 | Complete (normal) | None |
| 1 | Absent left PCoA | Left PCoA |
| 2 | Absent bilateral PCoA | Both PCoA |
| 3 | Absent left A1 (ACA) | Left ACA1 |
| 4 | Absent left P1 (PCA) | Left PCA1 |
| 5 | Absent right PCoA + left P1 | Right PCoA + Left PCA1 |

### Cerebral Flow Panel (bottom)
- 6 small D3 line charts in 2x3 grid (left hemisphere top row, right bottom)
- Labels: "MCA Left", "ACA Left", "PCA Left", etc.
- Shared time axis with aortic pressure, Y-axis in mL/s
- Compare mode overlays traces

## Information & Context
- Title + subtitle referencing Heldt (2002) and Ursino & Giannessi (2010)
- Collapsible "About" panel: model description, 57 ODEs, baroreflex, clinical significance
- All axes labeled with units
- CoW variant dropdown uses clinical names

## Tech Stack
- React + D3.js, Vite bundler
- Static site — no server required
- Deploy target: static hosting (Vercel / GitHub Pages / subdomain)
