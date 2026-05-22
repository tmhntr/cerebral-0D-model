# Cerebral Hemodynamic Simulator

A lumped-parameter model of coupled systemic and cerebral blood flow, developed as part of my Masters thesis. The model solves 57 coupled ODEs describing the cardiovascular system, cerebral vasculature (including the Circle of Willis), baroreflex control, and autoregulation.

An interactive web dashboard visualizes the simulation results, comparing healthy (normal sinus rhythm) and atrial fibrillation conditions across six Circle of Willis anatomical variants.

## Live Demo

[View the interactive dashboard](#) *(deploy URL TBD)*

## The Model

The simulator combines three published models:

- **Systemic circulation** (Heldt 2002) — Heart chambers, major arteries, venous return, and pulmonary circulation represented as lumped-parameter compartments
- **Cerebral vasculature** (Ursino & Giannessi 2010) — Circle of Willis anatomy with six cortical territories (MCA, ACA, PCA, left and right), autoregulation, and CO2 reactivity
- **Baroreflex** (Lin et al. 2012) — Autonomic control of heart rate, ventricular contractility, vascular tone, and venous volume

Atrial fibrillation is modeled by eliminating atrial contraction and introducing beat-to-beat variability (Scarsoglio et al. 2014).

### Circle of Willis Variants

The model supports six anatomical configurations:

| Code | Variant | Clinical Significance |
|------|---------|----------------------|
| 0 | Complete (normal) | All communicating arteries present |
| 1 | Absent left PCoA | Reduced posterior collateral on left |
| 2 | Absent bilateral PCoA | No posterior communication |
| 3 | Absent left A1 (ACA) | Left anterior territory fed via ACoA |
| 4 | Absent left P1 (PCA) | Left posterior territory fed via PCoA |
| 5 | Absent right PCoA + left P1 | Combined variant |

## Web Dashboard

The dashboard pre-computes 36 scenarios (2 conditions x 6 CoW variants x 3 heart rates) and displays:

- **Aortic pressure waveform** — Beat-to-beat pressure dynamics
- **Circle of Willis diagram** — SVG schematic with flow-encoded vessel coloring
- **Cerebral blood flow subplots** — Six territories showing regional perfusion differences

A compare mode overlays two scenarios for side-by-side analysis.

### Run Locally

```bash
cd web
npm install
npm run dev
```

Open http://localhost:5173.

### Regenerate Simulation Data

Requires SUNDIALS 7.x (`brew install sundials` on macOS):

```bash
make cbf
python3 scripts/precompute.py
```

This runs 36 simulations (~1 hour) and writes `web/public/data/scenarios.json`.

## C Simulator

### Dependencies

- C compiler (cc/gcc/clang)
- [SUNDIALS](https://computing.llnl.gov/projects/sundials) 7.x (ODE solver suite)
- Open MPI (linked by SUNDIALS)

On macOS: `brew install sundials`

### Build

```bash
make cbf
```

### Run

```bash
./cbf <run_index> <is_af> <cow_var> <hr_0>
```

Arguments:
- `run_index` — Simulation instance index (selects row from input files)
- `is_af` — 0 = normal sinus rhythm, 1 = atrial fibrillation
- `cow_var` — Circle of Willis variant (0-5)
- `hr_0` — Intrinsic heart rate (bpm)

Requires `input/` directory with `randomPars.dat`, `pinkNoise.dat`, and `expNoise.dat`.

### Deploy

```bash
docker build -t cerebral-sim .
docker run -p 8080:80 cerebral-sim
```

## License

GPL v3 — see [LICENSE](LICENSE).

## References

1. Heldt T. (2002). *Computational Models of Cardiovascular Response to Orthostatic Stress*. Journal of Applied Physiology.
2. Ursino M, Giannessi M. (2010). *A Model of Cerebrovascular Reactivity Including the Circle of Willis and Cortical Anastomoses*. Annals of Biomedical Engineering.
3. Lin J, et al. (2012). *A baroreflex model*. Proceedings of the Institution of Mechanical Engineers, Part H.
4. Scarsoglio S, et al. (2014). *Impact of atrial fibrillation on cerebral hemodynamics*.
