const COW_VARIANTS = [0, 1, 2, 3, 4, 5]
const HEART_RATES = [60, 75, 90]

function ScenarioSelector({ label, value, onChange, cowNames, className }) {
  return (
    <fieldset className={`scenario-selector ${className}`}>
      <legend>{label}</legend>
      <div className="selector-row">
        <div className="selector-group">
          <label>Condition</label>
          <select
            value={value.condition}
            onChange={(e) => onChange({ ...value, condition: e.target.value })}
          >
            <option value="nsr">NSR (Normal Sinus)</option>
            <option value="af">AF (Atrial Fibrillation)</option>
          </select>
        </div>

        <div className="selector-group">
          <label>Circle of Willis</label>
          <select
            value={value.cow}
            onChange={(e) => onChange({ ...value, cow: Number(e.target.value) })}
          >
            {COW_VARIANTS.map((v) => (
              <option key={v} value={v}>
                {v}: {cowNames[String(v)]}
              </option>
            ))}
          </select>
        </div>

        <div className="selector-group">
          <label>Heart Rate</label>
          <select
            value={value.hr}
            onChange={(e) => onChange({ ...value, hr: Number(e.target.value) })}
          >
            {HEART_RATES.map((hr) => (
              <option key={hr} value={hr}>
                {hr} bpm
              </option>
            ))}
          </select>
        </div>
      </div>
    </fieldset>
  )
}

export default function Controls({
  primary,
  setPrimary,
  compare,
  setCompare,
  compareEnabled,
  setCompareEnabled,
  cowNames,
}) {
  return (
    <div className="controls">
      <ScenarioSelector
        label="Primary Scenario"
        value={primary}
        onChange={setPrimary}
        cowNames={cowNames}
        className="primary"
      />

      <div className="compare-toggle">
        <input
          type="checkbox"
          id="compare-check"
          checked={compareEnabled}
          onChange={(e) => setCompareEnabled(e.target.checked)}
        />
        <label htmlFor="compare-check">Compare</label>
      </div>

      {compareEnabled && (
        <ScenarioSelector
          label="Compare Scenario"
          value={compare}
          onChange={setCompare}
          cowNames={cowNames}
          className="compare"
        />
      )}
    </div>
  )
}
