import { useState } from 'react'

export default function About() {
  const [open, setOpen] = useState(false)

  return (
    <>
      <button
        className="about-toggle"
        onClick={() => setOpen(!open)}
      >
        {open ? 'Hide' : 'About this model'}
      </button>

      {open && (
        <div className="about-panel">
          <p>
            This simulator models coupled systemic and cerebral hemodynamics
            using 57 ordinary differential equations.
          </p>
          <p>
            The systemic circulation is based on Heldt (2002), representing the
            heart, major arteries, and venous return with lumped-parameter
            compartments.
          </p>
          <p>
            The cerebral model follows Ursino &amp; Giannessi (2010),
            implementing the Circle of Willis, autoregulation, and CO
            <sub>2</sub> reactivity across six arterial territories (MCA, ACA,
            PCA, left and right).
          </p>
          <p>
            A baroreflex controller (Lin et al., 2012) modulates heart rate,
            ventricular contractility, and vascular tone.
          </p>
          <p>
            Atrial fibrillation is modeled by eliminating atrial contraction and
            introducing beat-to-beat variability following Scarsoglio et al.
            (2014).
          </p>
          <p>
            <a href="#" target="_blank" rel="noopener noreferrer">
              View source on GitHub
            </a>
          </p>
        </div>
      )}
    </>
  )
}
