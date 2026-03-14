# Statistical and Dynamical Analysis of Homogeneous Isotropic Turbulence

This repository contains a robust suite of tools for analysing **Homogeneous Isotropic Turbulence (HIT)** data. Using high-resolution simulation slices (1024³), the project explores the multi-scale nature of turbulence through Eulerian statistics, spectral analysis, and Lagrangian particle dynamics.

## Project Highlights
*   **Physics-Driven Insights:** Validated Kolmogorov’s $-5/3$ energy spectrum and structure function scaling laws.
*   **Computational Efficiency:** Vectorized NumPy operations and optimized spline interpolations for tracking $10^4$ particles.
*   **Advanced Visualisation:** High-fidelity contour plots, spectral energy density maps, and flow topology (Q-R) diagrams.
*   **Statistical Depth:** Characterised intermittency through high-order moments (Skewness/Kurtosis) and Probability Density Functions (PDFs).

---

## Technical Skill Set
*   **Languages:** Python
*   **Libraries:** NumPy, SciPy (Interpolation & Stats), Matplotlib, FFT (Fast Fourier Transforms)
*   **Mathematical Methods:** Spectral Analysis, Lagrangian Tracking, Tensor Calculus (Velocity Gradient Tensor), Numerical Differentiation.

---

## Phase 1: Eulerian Statistics & Field Characterisation
In this phase, the focus was on characterising the flow's "snapshot" properties. By calculating the velocity gradient tensor and enstrophy, I identified regions of high-intense rotation and dissipation.

### Key Results:
*   **Reynolds Number ($Re_\lambda$):** $\approx 7376$
*   **Intermittency:** Computed enstrophy kurtosis ($\approx 4385$), indicating the presence of extreme, localised events far from a normal distribution.
*   **Flow Visualisation:**
<table>
  <tr>
    <td><img src="Normalized_Vorticity.png" width="400" alt="Vorticity Field"/><br/><b>Normalized Vorticity Field</b></td>
    <td><img src="PDF_CDF_du_dx.png" width="400" alt="PDF du/dx"/><br/><b>Non-Gaussian Tails in Velocity Gradients</b></td>
  </tr>
</table>

---

## Phase 2: Spectral Analysis & Scaling Laws
Turbulence is defined by the transfer of energy across scales. This phase validates the **Kolmogorov-K41 theory**.

### Impactful Results:
*   **Energy Cascade:** Validated the $-5/3$ power law in both 1D and 2D spectra, identifying the inertial sub-range.
*   **Anomalous Scaling:** Measured the scaling exponents ($\tau_p$) for structure functions up to the 7th order. The deviation from the linear $p/3$ theory clearly demonstrates **internal intermittency**.

<p align="center">
  <img src="Structure Function Scaling Exponents.png" width="600" />
  <br/><i>Figure: Measured vs. Theoretical Scaling Exponents (highlighting intermittency effects).</i>
</p>

---

## Phase 3: Lagrangian Dynamics & Flow Topology
Using a custom-built particle tracker, I analysed how particles disperse in a chaotic flow field and characterised the local flow geometry.

### Key Results:
*   **Flow Topology (Q-R Plane):** Characterised the flow into "teardrop" shapes in the Q-R scatter plot, distinguishing between stable/unstable focal regions and strain-dominated zones.
*   **Richardson Pair Dispersion:** Analyzed particle separation. The growth rate $m \approx 3-4$ (Richardson's law) and positive **Lyapunov Exponents** ($\Lambda \approx 0.025$) provide a quantitative measure of the flow's chaotic sensitivity.

<table>
  <tr>
    <td><img src="Q_R_scatter.png" width="400" alt="QR Scatter"/><br/><b>Flow Topology (Q-R Scatter)</b></td>
    <td><img src="pair separation trajectories (Log-Log Scale).png" width="400" alt="Pair Separation"/><br/><b>Richardson Pair Dispersion Scaling</b></td>
  </tr>
</table>

---

## Numerical Summary
| Parameter | Value |
| :--- | :--- |
| **Reynolds Number (Re)** | 7376.65 |
| **Turbulent Diffusivity** | 1.172 |
| **Enstrophy Kurtosis** | 4385.6 |
| **Lyapunov Exponent ($\Lambda$)** | ~0.026 |
| **Grid resolution** | $1024^2$ slice |

---

## Repository Structure
*   `Turbulence_data_analysis-1.py`: Eulerian field analysis and basic statistics.
*   `Turbulence_data_analysis-2.py`: Spectral energy, correlations, and structure functions.
*   `Turbulence_data_analysis-3.py`: Particle tracking, pair dispersion, and Q-R topology.
