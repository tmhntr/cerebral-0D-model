import { useRef, useEffect } from 'react'
import * as d3 from 'd3'

const MARGIN = { top: 24, right: 12, bottom: 32, left: 40 }
const PRIMARY_COLOR = '#2563eb'
const COMPARE_COLOR = '#dc2626'

const FLOW_CHANNELS = [
  { col: 2, label: 'MCA Left' },
  { col: 3, label: 'ACA Left' },
  { col: 4, label: 'PCA Left' },
  { col: 5, label: 'MCA Right' },
  { col: 6, label: 'ACA Right' },
  { col: 7, label: 'PCA Right' },
]

function FlowSubplot({ primary, compare, channel }) {
  const svgRef = useRef()

  useEffect(() => {
    if (!primary) return

    const svg = d3.select(svgRef.current)
    svg.selectAll('*').remove()

    const container = svgRef.current.parentElement
    const width = container.clientWidth - 32
    const height = Math.min(width * 0.65, 200)

    svg.attr('viewBox', `0 0 ${width} ${height}`)

    const innerW = width - MARGIN.left - MARGIN.right
    const innerH = height - MARGIN.top - MARGIN.bottom

    const g = svg
      .append('g')
      .attr('transform', `translate(${MARGIN.left},${MARGIN.top})`)

    const primaryData = primary.data
    const compareData = compare ? compare.data : null
    const col = channel.col

    // Scales
    const xExtent = d3.extent(primaryData, (d) => d[0])
    if (compareData) {
      const cx = d3.extent(compareData, (d) => d[0])
      xExtent[0] = Math.min(xExtent[0], cx[0])
      xExtent[1] = Math.max(xExtent[1], cx[1])
    }

    let yMin = d3.min(primaryData, (d) => d[col])
    let yMax = d3.max(primaryData, (d) => d[col])
    if (compareData) {
      yMin = Math.min(yMin, d3.min(compareData, (d) => d[col]))
      yMax = Math.max(yMax, d3.max(compareData, (d) => d[col]))
    }
    const yPad = (yMax - yMin) * 0.15 || 0.5
    yMin -= yPad
    yMax += yPad

    const x = d3.scaleLinear().domain(xExtent).range([0, innerW])
    const y = d3.scaleLinear().domain([yMin, yMax]).range([innerH, 0])

    // Axes
    g.append('g')
      .attr('class', 'axis')
      .attr('transform', `translate(0,${innerH})`)
      .call(d3.axisBottom(x).ticks(4))

    g.append('g')
      .attr('class', 'axis')
      .call(d3.axisLeft(y).ticks(4))

    // Axis labels
    g.append('text')
      .attr('class', 'axis-label')
      .attr('x', innerW / 2)
      .attr('y', innerH + 28)
      .attr('text-anchor', 'middle')
      .text('Time (s)')

    g.append('text')
      .attr('class', 'axis-label')
      .attr('transform', 'rotate(-90)')
      .attr('x', -innerH / 2)
      .attr('y', -30)
      .attr('text-anchor', 'middle')
      .text('mL/s')

    // Title
    svg
      .append('text')
      .attr('class', 'chart-title')
      .attr('x', width / 2)
      .attr('y', 14)
      .attr('text-anchor', 'middle')
      .text(channel.label)

    // Line generator
    const line = d3
      .line()
      .x((d) => x(d[0]))
      .y((d) => y(d[col]))
      .curve(d3.curveMonotoneX)

    // Primary trace
    g.append('path')
      .datum(primaryData)
      .attr('fill', 'none')
      .attr('stroke', PRIMARY_COLOR)
      .attr('stroke-width', 1.5)
      .attr('d', line)

    // Compare trace
    if (compareData) {
      g.append('path')
        .datum(compareData)
        .attr('fill', 'none')
        .attr('stroke', COMPARE_COLOR)
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '6,3')
        .attr('d', line)
    }
  }, [primary, compare, channel])

  return (
    <div className="chart-panel">
      <svg ref={svgRef} preserveAspectRatio="xMidYMid meet" />
    </div>
  )
}

export default function FlowCharts({ primary, compare }) {
  return (
    <div className="flow-grid">
      {FLOW_CHANNELS.map((ch) => (
        <FlowSubplot
          key={ch.label}
          primary={primary}
          compare={compare}
          channel={ch}
        />
      ))}
    </div>
  )
}
