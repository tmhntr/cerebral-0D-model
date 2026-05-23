import { useState, useEffect } from 'react'
import Controls from './components/Controls'
import PressureChart from './components/PressureChart'
import FlowCharts from './components/FlowCharts'
import CowDiagram from './components/CowDiagram'
import About from './components/About'
import './App.css'

const DEFAULT_STATE = { condition: 'nsr', cow: 0, hr: 75 }

function buildKey({ condition, cow, hr }) {
  return `${condition}_cow${cow}_hr${hr}`
}

export default function App() {
  const [scenarioData, setScenarioData] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  const [primary, setPrimary] = useState(DEFAULT_STATE)
  const [compare, setCompare] = useState(DEFAULT_STATE)
  const [compareEnabled, setCompareEnabled] = useState(false)

  useEffect(() => {
    fetch(`${import.meta.env.BASE_URL}data/scenarios.json`)
      .then((res) => {
        if (!res.ok) throw new Error(`HTTP ${res.status}`)
        return res.json()
      })
      .then((data) => {
        setScenarioData(data)
        setLoading(false)
      })
      .catch((err) => {
        setError(err.message)
        setLoading(false)
      })
  }, [])

  if (loading) return <div className="loading">Loading simulation data...</div>
  if (error) return <div className="error">Failed to load data: {error}</div>

  const primaryKey = buildKey(primary)
  const compareKey = buildKey(compare)
  const primaryScenario = scenarioData.scenarios[primaryKey]
  const compareScenario = compareEnabled
    ? scenarioData.scenarios[compareKey]
    : null

  return (
    <div className="app">
      <header className="header">
        <h1>Cerebral Blood Flow Simulator</h1>
        <p className="subtitle">
          Coupled systemic-cerebral hemodynamics with Circle of Willis variants
        </p>
        <About />
      </header>

      <Controls
        primary={primary}
        setPrimary={setPrimary}
        compare={compare}
        setCompare={setCompare}
        compareEnabled={compareEnabled}
        setCompareEnabled={setCompareEnabled}
        cowNames={scenarioData.cow_names}
      />

      {!primaryScenario ? (
        <div className="error">
          Scenario not found: {primaryKey}
        </div>
      ) : (
        <>
          <div className="top-row">
            <PressureChart
              primary={primaryScenario}
              compare={compareScenario}
            />
            <CowDiagram
              primary={primaryScenario}
              compare={compareScenario}
              cowVariant={primary.cow}
            />
          </div>

          <FlowCharts
            primary={primaryScenario}
            compare={compareScenario}
          />
        </>
      )}
    </div>
  )
}
