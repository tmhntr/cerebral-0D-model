# Web Visualization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build an interactive static web dashboard that visualizes pre-computed cerebral blood flow simulation results, with Circle of Willis diagram, pressure waveforms, and flow subplots.

**Architecture:** Python script runs the C simulator across 36 scenario combinations and extracts 10 steady-state beats per scenario into a JSON file. React + D3 frontend loads this JSON and renders three panels (pressure waveform, CoW SVG diagram, 6 flow subplots) with instant scenario switching and compare mode.

**Tech Stack:** Python 3 (precompute), React 18, D3.js v7, Vite 5

---

### Task 1: Precompute Script

**Files:**
- Create: `scripts/precompute.py`

**Step 1: Write the precompute script**

This script runs `./cbf` for each of 36 scenarios, parses the output, extracts beats 100-109, and writes JSON.

```python
#!/usr/bin/env python3
"""Pre-compute simulation scenarios for the web visualization."""

import subprocess
import json
import os
import sys

SCENARIOS = {
    "conditions": [0, 1],  # 0=NSR, 1=AF
    "cow_variants": [0, 1, 2, 3, 4, 5],
    "heart_rates": [60, 75, 90],
}

CONDITION_NAMES = {0: "nsr", 1: "af"}
COW_NAMES = {
    0: "Complete (normal)",
    1: "Absent left PCoA",
    2: "Absent bilateral PCoA",
    3: "Absent left A1 (ACA)",
    4: "Absent left P1 (PCA)",
    5: "Absent right PCoA + left P1",
}
COLUMNS = ["time", "P_a", "q_ml", "q_al", "q_pl", "q_mr", "q_ar", "q_pr"]
BEAT_START = 100  # first beat to extract (after transients settle)
BEAT_END = 110    # exclusive

CBF_BINARY = os.path.join(os.path.dirname(__file__), "..", "cbf")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "web", "public", "data")


def run_scenario(run_index, is_af, cow_var, hr):
    """Run the C simulator and return (states_file, beat_times_file)."""
    subprocess.run(
        [CBF_BINARY, str(run_index), str(is_af), str(cow_var), str(hr)],
        check=True,
        cwd=os.path.join(os.path.dirname(__file__), ".."),
    )
    cond_label = "AF" if is_af else "NSR"
    states_file = f"statesOutput.{run_index:05d}.{cond_label}.dat"
    beats_file = f"endDiastolicTime.{run_index:05d}.{cond_label}.dat"
    return states_file, beats_file


def extract_beats(states_path, beats_path, beat_start, beat_end):
    """Extract timesteps between beat_start and beat_end."""
    with open(beats_path) as f:
        beat_times = [float(line.strip()) for line in f if line.strip()]

    t_start = beat_times[beat_start - 1]  # beat_start is 1-indexed in file
    t_end = beat_times[beat_end - 1]

    rows = []
    with open(states_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            t = float(parts[0])
            if t < t_start:
                continue
            if t >= t_end:
                break
            # Normalize time to start at 0
            row = [round(t - t_start, 4)] + [round(float(v), 6) for v in parts[1:]]
            rows.append(row)
    return rows


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    base_dir = os.path.join(os.path.dirname(__file__), "..")

    result = {"columns": COLUMNS, "cow_names": COW_NAMES, "scenarios": {}}
    run_index = 0

    total = (len(SCENARIOS["conditions"])
             * len(SCENARIOS["cow_variants"])
             * len(SCENARIOS["heart_rates"]))
    done = 0

    for is_af in SCENARIOS["conditions"]:
        for cow in SCENARIOS["cow_variants"]:
            for hr in SCENARIOS["heart_rates"]:
                key = f"{CONDITION_NAMES[is_af]}_cow{cow}_hr{hr}"
                done += 1
                print(f"[{done}/{total}] Running {key}...")

                states_file, beats_file = run_scenario(run_index, is_af, cow, hr)
                states_path = os.path.join(base_dir, states_file)
                beats_path = os.path.join(base_dir, beats_file)

                rows = extract_beats(states_path, beats_path, BEAT_START, BEAT_END)
                result["scenarios"][key] = {
                    "condition": CONDITION_NAMES[is_af],
                    "cow": cow,
                    "hr": hr,
                    "data": rows,
                }

                # Clean up output files
                os.remove(states_path)
                os.remove(beats_path)
                params_file = os.path.join(base_dir, f"parameters.{run_index:05d}.dat")
                if os.path.exists(params_file):
                    os.remove(params_file)

                run_index += 1

    out_path = os.path.join(OUTPUT_DIR, "scenarios.json")
    with open(out_path, "w") as f:
        json.dump(result, f)

    size_mb = os.path.getsize(out_path) / (1024 * 1024)
    print(f"Wrote {out_path} ({size_mb:.1f} MB, {len(result['scenarios'])} scenarios)")


if __name__ == "__main__":
    main()
```

**Step 2: Run it (takes ~30-60 min for 36 scenarios)**

```bash
cd /Users/timhunter/Developer/revival/repos/cerebral-0D-model
python3 scripts/precompute.py
```

Expected: `web/public/data/scenarios.json` created with 36 scenario keys.

**Step 3: Commit**

```bash
git add scripts/precompute.py
git commit -m "feat: add precompute script for 36 simulation scenarios"
```

---

### Task 2: Scaffold React + Vite App

**Files:**
- Create: `web/` directory (via `npm create vite`)
- Create: `web/src/App.jsx`
- Create: `web/src/main.jsx`
- Create: `web/src/index.css`

**Step 1: Scaffold and install deps**

```bash
cd /Users/timhunter/Developer/revival/repos/cerebral-0D-model
npm create vite@latest web -- --template react
cd web
npm install
npm install d3
```

**Step 2: Clean up scaffolded files**

Remove default Vite boilerplate from `App.jsx`, `App.css`, etc. Replace with minimal shell:

`web/src/App.jsx`:
```jsx
import { useState, useEffect } from 'react'
import './App.css'

const COW_NAMES = {
  0: "Complete (normal)",
  1: "Absent left PCoA",
  2: "Absent bilateral PCoA",
  3: "Absent left A1 (ACA)",
  4: "Absent left P1 (PCA)",
  5: "Absent right PCoA + left P1",
}

function App() {
  const [data, setData] = useState(null)
  const [condition, setCondition] = useState('nsr')
  const [cow, setCow] = useState(0)
  const [hr, setHr] = useState(75)
  const [compare, setCompare] = useState(false)
  const [condition2, setCondition2] = useState('af')
  const [cow2, setCow2] = useState(0)
  const [hr2, setHr2] = useState(75)

  useEffect(() => {
    fetch('/data/scenarios.json')
      .then(r => r.json())
      .then(setData)
  }, [])

  if (!data) return <div className="loading">Loading simulation data...</div>

  const key1 = `${condition}_cow${cow}_hr${hr}`
  const key2 = `${condition2}_cow${cow2}_hr${hr2}`
  const scenario1 = data.scenarios[key1]
  const scenario2 = compare ? data.scenarios[key2] : null

  return (
    <div className="app">
      <header>
        <h1>Cerebral Hemodynamic Simulator</h1>
        <p className="subtitle">
          Lumped-parameter model of systemic and cerebral blood flow
          — Heldt (2002) &amp; Ursino &amp; Giannessi (2010)
        </p>
      </header>
      <Controls
        condition={condition} setCondition={setCondition}
        cow={cow} setCow={setCow}
        hr={hr} setHr={setHr}
        compare={compare} setCompare={setCompare}
        condition2={condition2} setCondition2={setCondition2}
        cow2={cow2} setCow2={setCow2}
        hr2={hr2} setHr2={setHr2}
      />
      <div className="panels">
        <div className="top-row">
          <PressureChart scenario1={scenario1} scenario2={scenario2} columns={data.columns} />
          <CowDiagram scenario1={scenario1} scenario2={scenario2} cow={cow} cow2={compare ? cow2 : null} />
        </div>
        <FlowCharts scenario1={scenario1} scenario2={scenario2} columns={data.columns} />
      </div>
    </div>
  )
}
```

**Step 3: Verify dev server starts**

```bash
cd web && npm run dev
```

Expected: Vite dev server on http://localhost:5173, shows loading state.

**Step 4: Commit**

```bash
git add web/
git commit -m "feat: scaffold React + Vite app with app shell and controls"
```

---

### Task 3: Controls Component

**Files:**
- Create: `web/src/components/Controls.jsx`
- Create: `web/src/components/Controls.css`

**Step 1: Build the Controls component**

Two side-by-side scenario selectors. Left always active, right active only when compare is on. Each has: Condition dropdown (NSR/AF), CoW variant dropdown (clinical names), HR dropdown (60/75/90).

```jsx
import './Controls.css'

const COW_NAMES = {
  0: "Complete (normal)",
  1: "Absent left PCoA",
  2: "Absent bilateral PCoA",
  3: "Absent left A1 (ACA)",
  4: "Absent left P1 (PCA)",
  5: "Absent right PCoA + left P1",
}

export default function Controls({
  condition, setCondition, cow, setCow, hr, setHr,
  compare, setCompare,
  condition2, setCondition2, cow2, setCow2, hr2, setHr2,
}) {
  return (
    <div className="controls">
      <ScenarioSelector
        label="Scenario"
        condition={condition} setCondition={setCondition}
        cow={cow} setCow={setCow}
        hr={hr} setHr={setHr}
        color="#2563eb"
      />
      <label className="compare-toggle">
        <input type="checkbox" checked={compare} onChange={e => setCompare(e.target.checked)} />
        Compare
      </label>
      {compare && (
        <ScenarioSelector
          label="Compare with"
          condition={condition2} setCondition={setCondition2}
          cow={cow2} setCow={setCow2}
          hr={hr2} setHr={setHr2}
          color="#dc2626"
        />
      )}
    </div>
  )
}

function ScenarioSelector({ label, condition, setCondition, cow, setCow, hr, setHr, color }) {
  return (
    <fieldset className="scenario-selector" style={{ borderColor: color }}>
      <legend style={{ color }}>{label}</legend>
      <label>
        Condition
        <select value={condition} onChange={e => setCondition(e.target.value)}>
          <option value="nsr">Normal Sinus Rhythm</option>
          <option value="af">Atrial Fibrillation</option>
        </select>
      </label>
      <label>
        Circle of Willis
        <select value={cow} onChange={e => setCow(Number(e.target.value))}>
          {Object.entries(COW_NAMES).map(([k, v]) => (
            <option key={k} value={k}>{v}</option>
          ))}
        </select>
      </label>
      <label>
        Heart Rate
        <select value={hr} onChange={e => setHr(Number(e.target.value))}>
          {[60, 75, 90].map(h => (
            <option key={h} value={h}>{h} bpm</option>
          ))}
        </select>
      </label>
    </fieldset>
  )
}
```

**Step 2: Verify controls render and update state**

Open dev server, confirm dropdowns appear and switching them updates the page.

**Step 3: Commit**

```bash
git add web/src/components/
git commit -m "feat: add scenario controls with compare mode"
```

---

### Task 4: Aortic Pressure Waveform Chart

**Files:**
- Create: `web/src/components/PressureChart.jsx`

**Step 1: Build the D3 line chart component**

```jsx
import { useRef, useEffect } from 'react'
import * as d3 from 'd3'

const MARGIN = { top: 20, right: 20, bottom: 40, left: 55 }

export default function PressureChart({ scenario1, scenario2, columns }) {
  const svgRef = useRef()

  useEffect(() => {
    if (!scenario1) return
    const svg = d3.select(svgRef.current)
    svg.selectAll('*').remove()

    const width = svgRef.current.clientWidth
    const height = svgRef.current.clientHeight
    const w = width - MARGIN.left - MARGIN.right
    const h = height - MARGIN.top - MARGIN.bottom

    const g = svg.append('g').attr('transform', `translate(${MARGIN.left},${MARGIN.top})`)

    const tIdx = 0 // time column
    const pIdx = 1 // P_a column

    const allData = scenario1.data.concat(scenario2 ? scenario2.data : [])
    const xExtent = d3.extent(scenario1.data, d => d[tIdx])
    const yExtent = d3.extent(allData, d => d[pIdx])
    yExtent[0] = Math.min(yExtent[0], 30)
    yExtent[1] = Math.max(yExtent[1], 130)

    const x = d3.scaleLinear().domain(xExtent).range([0, w])
    const y = d3.scaleLinear().domain(yExtent).nice().range([h, 0])

    // Axes
    g.append('g').attr('transform', `translate(0,${h})`).call(d3.axisBottom(x).ticks(6))
    g.append('g').call(d3.axisLeft(y).ticks(6))

    // Axis labels
    g.append('text').attr('x', w / 2).attr('y', h + 35).attr('text-anchor', 'middle')
      .attr('class', 'axis-label').text('Time (s)')
    g.append('text').attr('transform', 'rotate(-90)').attr('y', -45).attr('x', -h / 2)
      .attr('text-anchor', 'middle').attr('class', 'axis-label').text('Pressure (mmHg)')

    // Title
    g.append('text').attr('x', w / 2).attr('y', -5).attr('text-anchor', 'middle')
      .attr('class', 'chart-title').text('Aortic Pressure')

    const line = d3.line().x(d => x(d[tIdx])).y(d => y(d[pIdx]))

    // Primary trace
    g.append('path').datum(scenario1.data)
      .attr('fill', 'none').attr('stroke', '#2563eb').attr('stroke-width', 1.5)
      .attr('d', line)

    // Compare trace
    if (scenario2) {
      g.append('path').datum(scenario2.data)
        .attr('fill', 'none').attr('stroke', '#dc2626').attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '6,3').attr('d', line)
    }
  }, [scenario1, scenario2])

  return <svg ref={svgRef} className="pressure-chart" />
}
```

**Step 2: Verify chart renders with real data**

Load the dev server with `scenarios.json` in `web/public/data/`. Confirm pressure waveform appears.

**Step 3: Commit**

```bash
git add web/src/components/PressureChart.jsx
git commit -m "feat: add aortic pressure waveform chart"
```

---

### Task 5: Cerebral Flow Subplots

**Files:**
- Create: `web/src/components/FlowCharts.jsx`

**Step 1: Build the 2x3 flow chart grid**

Similar to PressureChart but renders 6 small charts. Column indices:
- q_ml=2, q_al=3, q_pl=4 (left hemisphere)
- q_mr=5, q_ar=6, q_pr=7 (right hemisphere)

```jsx
import { useRef, useEffect } from 'react'
import * as d3 from 'd3'

const FLOWS = [
  { idx: 2, label: 'MCA Left' },
  { idx: 3, label: 'ACA Left' },
  { idx: 4, label: 'PCA Left' },
  { idx: 5, label: 'MCA Right' },
  { idx: 6, label: 'ACA Right' },
  { idx: 7, label: 'PCA Right' },
]

const MARGIN = { top: 25, right: 10, bottom: 30, left: 45 }

export default function FlowCharts({ scenario1, scenario2, columns }) {
  const containerRef = useRef()

  useEffect(() => {
    if (!scenario1 || !containerRef.current) return
    const container = d3.select(containerRef.current)
    container.selectAll('*').remove()

    FLOWS.forEach((flow, i) => {
      const svg = container.append('svg').attr('class', 'flow-subplot')
      const width = svg.node().clientWidth
      const height = svg.node().clientHeight
      const w = width - MARGIN.left - MARGIN.right
      const h = height - MARGIN.top - MARGIN.bottom

      const g = svg.append('g').attr('transform', `translate(${MARGIN.left},${MARGIN.top})`)

      const tIdx = 0
      const xExtent = d3.extent(scenario1.data, d => d[tIdx])
      const allData = scenario1.data.concat(scenario2 ? scenario2.data : [])
      const yExtent = d3.extent(allData, d => d[flow.idx])

      const x = d3.scaleLinear().domain(xExtent).range([0, w])
      const y = d3.scaleLinear().domain(yExtent).nice().range([h, 0])

      g.append('g').attr('transform', `translate(0,${h})`).call(d3.axisBottom(x).ticks(4))
      g.append('g').call(d3.axisLeft(y).ticks(4))

      // Y-axis label only on leftmost column
      if (i % 3 === 0) {
        g.append('text').attr('transform', 'rotate(-90)').attr('y', -35).attr('x', -h / 2)
          .attr('text-anchor', 'middle').attr('class', 'axis-label').text('mL/s')
      }

      g.append('text').attr('x', w / 2).attr('y', -8).attr('text-anchor', 'middle')
        .attr('class', 'chart-title').text(flow.label)

      const line = d3.line().x(d => x(d[tIdx])).y(d => y(d[flow.idx]))

      g.append('path').datum(scenario1.data)
        .attr('fill', 'none').attr('stroke', '#2563eb').attr('stroke-width', 1.5).attr('d', line)

      if (scenario2) {
        g.append('path').datum(scenario2.data)
          .attr('fill', 'none').attr('stroke', '#dc2626').attr('stroke-width', 1.5)
          .attr('stroke-dasharray', '6,3').attr('d', line)
      }
    })
  }, [scenario1, scenario2])

  return <div ref={containerRef} className="flow-charts" />
}
```

**Step 2: Verify 6 subplots render**

**Step 3: Commit**

```bash
git add web/src/components/FlowCharts.jsx
git commit -m "feat: add 6 cerebral flow subplot charts"
```

---

### Task 6: Circle of Willis SVG Diagram

**Files:**
- Create: `web/src/components/CowDiagram.jsx`
- Create: `web/src/components/CowDiagram.css`

**Step 1: Build the anatomical SVG diagram**

This is the centerpiece. Hand-drawn SVG paths for the Circle of Willis vessels, positioned in a top-down anatomical layout. Vessels are colored by mean flow and missing vessels are shown as dashed grey.

The component needs to:
1. Define SVG path coordinates for each vessel segment
2. Compute mean flows from scenario data
3. Map flow magnitude to stroke width + color
4. Apply dashed-grey style for missing vessels based on CoW variant
5. Add hover tooltips

Key vessel paths (approximate coordinates in a 300x350 SVG viewBox):
- Left/Right ICA: vertical lines rising from bottom-left/right
- BA: vertical line rising from bottom-center
- Left/Right MCA: branches going left/right from ICA tops
- ACA1 L/R: branches going up-inward from ICA tops
- ACoA: horizontal connector between ACA1 L/R
- ACA2 L/R: continuing upward from ACoA junction
- PCA1 L/R: branches from BA top going left/right
- PCA2 L/R: continuing outward from PCA1
- PCoA L/R: connectors between ICA and PCA1

**Step 2: Define which vessels are absent per CoW variant**

```javascript
const MISSING_VESSELS = {
  0: [],
  1: ['PCoA_L'],
  2: ['PCoA_L', 'PCoA_R'],
  3: ['ACA1_L'],
  4: ['PCA1_L'],
  5: ['PCoA_R', 'PCA1_L'],
}
```

**Step 3: Compute mean flows from scenario data for coloring**

The scenario data has q_ml, q_al, q_pl, q_mr, q_ar, q_pr (column indices 2-7). These are the distal cerebral flows. For the CoW vessels (ICA, BA, PCoA, ACoA, MCA, ACA, PCA), we approximate:
- ICA_L flow ≈ q_ml + q_al (feeds left MCA + ACA territories)
- ICA_R flow ≈ q_mr + q_ar
- BA flow ≈ q_pl + q_pr (feeds posterior territories)
- MCA_L/R = q_ml / q_mr
- ACA_L/R = q_al / q_ar
- PCA_L/R = q_pl / q_pr

**Step 4: Verify diagram renders with correct anatomy and coloring**

**Step 5: Commit**

```bash
git add web/src/components/CowDiagram.*
git commit -m "feat: add Circle of Willis SVG diagram with flow coloring"
```

---

### Task 7: Styling and Layout

**Files:**
- Modify: `web/src/App.css`
- Modify: `web/src/index.css`

**Step 1: Apply CSS grid layout**

```css
.app {
  max-width: 1200px;
  margin: 0 auto;
  padding: 1rem;
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  color: #1a1a2e;
  background: #f8f9fa;
}

header { text-align: center; margin-bottom: 1rem; }
header h1 { margin: 0; font-size: 1.6rem; }
.subtitle { color: #666; font-size: 0.85rem; margin: 0.25rem 0 0; }

.top-row {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 1rem;
  margin-bottom: 1rem;
}

.pressure-chart { width: 100%; height: 300px; }

.flow-charts {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 0.5rem;
}
.flow-subplot { width: 100%; height: 200px; }
```

**Step 2: Verify responsive layout**

**Step 3: Commit**

```bash
git add web/src/App.css web/src/index.css
git commit -m "feat: add dashboard layout and styling"
```

---

### Task 8: About Panel

**Files:**
- Create: `web/src/components/About.jsx`

**Step 1: Add collapsible info panel**

Toggle button in header. Content: model description, 57 ODEs, references.

**Step 2: Commit**

```bash
git add web/src/components/About.jsx
git commit -m "feat: add collapsible About panel with model description"
```

---

### Task 9: Run Precompute and Verify End-to-End

**Step 1: Build the C simulator**

```bash
cd /Users/timhunter/Developer/revival/repos/cerebral-0D-model
make clean && make cbf
```

**Step 2: Run precompute (takes ~30-60 min)**

```bash
python3 scripts/precompute.py
```

**Step 3: Verify JSON output**

```bash
python3 -c "import json; d=json.load(open('web/public/data/scenarios.json')); print(len(d['scenarios']), 'scenarios'); print(list(d['scenarios'].keys())[:5])"
```

Expected: `36 scenarios` and keys like `nsr_cow0_hr60`, `nsr_cow0_hr75`, etc.

**Step 4: Start dev server and verify full dashboard**

```bash
cd web && npm run dev
```

Open http://localhost:5173, verify:
- Pressure waveform renders with oscillating pressure
- 6 flow subplots show cerebral blood flows
- CoW diagram shows anatomy with flow colors
- Switching dropdowns instantly updates all panels
- Compare mode overlays two traces

**Step 5: Commit data file**

```bash
git add web/public/data/scenarios.json
git commit -m "feat: add pre-computed scenario data (36 scenarios)"
```

---

### Task 10: Production Build and Cleanup

**Step 1: Build for production**

```bash
cd web && npm run build
```

Verify `web/dist/` contains the static site.

**Step 2: Test production build**

```bash
npx serve web/dist
```

Open in browser, verify all features work.

**Step 3: Final commit**

```bash
git add -A
git commit -m "feat: complete web visualization dashboard"
```
