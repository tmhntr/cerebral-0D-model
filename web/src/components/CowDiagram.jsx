import { useRef, useEffect, useState } from 'react'
import * as d3 from 'd3'

const WIDTH = 400
const HEIGHT = 400
const CX = WIDTH / 2
const CY = HEIGHT / 2

// Missing vessels per CoW variant
const MISSING = {
  0: [],
  1: ['PCoA_L'],
  2: ['PCoA_L', 'PCoA_R'],
  3: ['ACA1_L'],
  4: ['PCA1_L'],
  5: ['PCoA_R', 'PCA1_L'],
}

// Vessel path definitions: [x1, y1, x2, y2]
const VESSELS = {
  // Basilar artery - comes up from bottom center
  BA:       { path: [CX, CY + 120, CX, CY + 40],   label: [CX + 14, CY + 85] },

  // BA splits into left and right PCA1
  PCA1_L:   { path: [CX, CY + 40,  CX - 55, CY + 20],  label: [CX - 45, CY + 15] },
  PCA1_R:   { path: [CX, CY + 40,  CX + 55, CY + 20],  label: [CX + 25, CY + 15] },

  // PCA2 continues laterally/up
  PCA2_L:   { path: [CX - 55, CY + 20, CX - 110, CY - 10], label: [CX - 120, CY + 8] },
  PCA2_R:   { path: [CX + 55, CY + 20, CX + 110, CY - 10], label: [CX + 80, CY + 8] },

  // PCoA connects PCA1 junction to ICA
  PCoA_L:   { path: [CX - 55, CY + 20, CX - 55, CY - 30],  label: [CX - 78, CY - 5] },
  PCoA_R:   { path: [CX + 55, CY + 20, CX + 55, CY - 30],  label: [CX + 60, CY - 5] },

  // ICAs come up from bottom-left/right to PCoA junction level
  ICA_L:    { path: [CX - 80, CY + 120, CX - 55, CY - 30],  label: [CX - 90, CY + 50] },
  ICA_R:    { path: [CX + 80, CY + 120, CX + 55, CY - 30],  label: [CX + 65, CY + 50] },

  // MCA goes laterally outward from ICA top
  MCA_L:    { path: [CX - 55, CY - 30, CX - 130, CY - 60],  label: [CX - 140, CY - 55] },
  MCA_R:    { path: [CX + 55, CY - 30, CX + 130, CY - 60],  label: [CX + 100, CY - 55] },

  // ACA1 goes medially upward from ICA top
  ACA1_L:   { path: [CX - 55, CY - 30, CX - 20, CY - 65],   label: [CX - 50, CY - 55] },
  ACA1_R:   { path: [CX + 55, CY - 30, CX + 20, CY - 65],   label: [CX + 25, CY - 55] },

  // ACoA connects left ACA1 to right ACA1 (horizontal)
  ACoA:     { path: [CX - 20, CY - 65, CX + 20, CY - 65],    label: [CX - 8, CY - 72] },

  // ACA2 continues upward from ACoA junction
  ACA2_L:   { path: [CX - 20, CY - 65, CX - 25, CY - 130],  label: [CX - 50, CY - 105] },
  ACA2_R:   { path: [CX + 20, CY - 65, CX + 25, CY - 130],  label: [CX + 30, CY - 105] },
}

// Map computed flow values to vessels for stroke encoding
function computeFlows(data) {
  if (!data || data.length === 0) return {}
  const mean = (col) => d3.mean(data, (d) => d[col])
  const mca_l = mean(2)
  const aca_l = mean(3)
  const pca_l = mean(4)
  const mca_r = mean(5)
  const aca_r = mean(6)
  const pca_r = mean(7)
  return {
    BA: pca_l + pca_r,
    PCA1_L: pca_l,
    PCA1_R: pca_r,
    PCA2_L: pca_l,
    PCA2_R: pca_r,
    PCoA_L: 0.1,   // communicating arteries carry small flow normally
    PCoA_R: 0.1,
    ICA_L: mca_l + aca_l,
    ICA_R: mca_r + aca_r,
    MCA_L: mca_l,
    MCA_R: mca_r,
    ACA1_L: aca_l,
    ACA1_R: aca_r,
    ACoA: 0.05,
    ACA2_L: aca_l,
    ACA2_R: aca_r,
  }
}

export default function CowDiagram({ primary, compare, cowVariant }) {
  const svgRef = useRef()
  const [tooltip, setTooltip] = useState(null)

  useEffect(() => {
    if (!primary) return

    const svg = d3.select(svgRef.current)
    svg.selectAll('*').remove()

    svg.attr('viewBox', `0 0 ${WIDTH} ${HEIGHT}`)

    const flows = computeFlows(primary.data)
    const compareFlows = compare ? computeFlows(compare.data) : null
    const missing = MISSING[cowVariant] || []

    // Collect all flow values for scaling
    const allFlows = Object.values(flows).filter((v) => v > 0)
    const flowMin = d3.min(allFlows) || 0.01
    const flowMax = d3.max(allFlows) || 5

    const strokeScale = d3
      .scaleLinear()
      .domain([flowMin, flowMax])
      .range([2, 8])
      .clamp(true)

    // Color scale: low flow = blue, high flow = red
    const colorScale = d3
      .scaleSequential(d3.interpolateRdYlBu)
      .domain([flowMax, flowMin]) // reversed so high = red

    // Title
    svg
      .append('text')
      .attr('class', 'chart-title')
      .attr('x', CX)
      .attr('y', 16)
      .attr('text-anchor', 'middle')
      .text('Circle of Willis')

    const g = svg.append('g').attr('transform', `translate(0, 20)`)

    // Draw vessels
    Object.entries(VESSELS).forEach(([name, vessel]) => {
      const isMissing = missing.includes(name)
      const flow = flows[name] || 0
      const [x1, y1, x2, y2] = vessel.path

      // Vessel line
      g.append('line')
        .attr('x1', x1)
        .attr('y1', y1)
        .attr('x2', x2)
        .attr('y2', y2)
        .attr('stroke', isMissing ? '#ccc' : colorScale(flow))
        .attr('stroke-width', isMissing ? 1.5 : strokeScale(flow))
        .attr('stroke-dasharray', isMissing ? '4,3' : 'none')
        .attr('stroke-linecap', 'round')
        .attr('cursor', 'pointer')
        .on('mouseenter', function (event) {
          const rect = svgRef.current.getBoundingClientRect()
          const svgX = event.clientX - rect.left
          const svgY = event.clientY - rect.top
          setTooltip({
            x: svgX,
            y: svgY - 30,
            text: isMissing
              ? `${name} (absent)`
              : `${name}: ${flow.toFixed(2)} mL/s`,
          })
          d3.select(this).attr('stroke-width', isMissing ? 3 : strokeScale(flow) + 2)
        })
        .on('mouseleave', function () {
          setTooltip(null)
          d3.select(this).attr('stroke-width', isMissing ? 1.5 : strokeScale(flow))
        })

      // Label
      const [lx, ly] = vessel.label
      // Simplify label names for display
      const displayName = name
        .replace('_L', ' L')
        .replace('_R', ' R')
        .replace(/[12]$/, (m) => '\u2081\u2082'[m - 1] || m)

      g.append('text')
        .attr('x', lx)
        .attr('y', ly)
        .attr('font-size', 9)
        .attr('fill', isMissing ? '#bbb' : '#555')
        .attr('text-anchor', 'start')
        .text(displayName)
    })

    // Draw junction circles at key branch points
    const junctions = [
      [CX, CY + 40],       // BA bifurcation
      [CX - 55, CY + 20],  // Left PCA1/PCoA junction
      [CX + 55, CY + 20],  // Right PCA1/PCoA junction
      [CX - 55, CY - 30],  // Left ICA top (MCA/ACA branch)
      [CX + 55, CY - 30],  // Right ICA top
      [CX - 20, CY - 65],  // Left ACoA junction
      [CX + 20, CY - 65],  // Right ACoA junction
    ]

    junctions.forEach(([jx, jy]) => {
      g.append('circle')
        .attr('cx', jx)
        .attr('cy', jy)
        .attr('r', 3)
        .attr('fill', '#374151')
    })
  }, [primary, compare, cowVariant])

  return (
    <div className="chart-panel cow-panel">
      <svg ref={svgRef} preserveAspectRatio="xMidYMid meet" />
      {tooltip && (
        <div
          className="cow-tooltip"
          style={{ left: tooltip.x, top: tooltip.y }}
        >
          {tooltip.text}
        </div>
      )}
    </div>
  )
}
