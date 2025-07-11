import math

def evaluate_quantum_correction():
    """
    Calculates the quantum correction to conductivity (weak localization)
    for an electron in a bulk (3D) semiconductor.
    """
    # --- 1. Define Physical Constants ---
    e = 1.60217663e-19  # Elementary charge (Coulombs)
    hbar = 1.05457182e-34 # Reduced Planck constant (Joule-seconds)
    m_e = 9.1093837e-31   # Electron rest mass (kg)
    
    # --- 2. Define Typical Semiconductor Parameters ---
    # Using parameters representative of a doped n-type GaAs sample at low temperature.
    # Effective mass ratio for GaAs
    m_star_ratio = 0.067
    # Carrier (electron) density in m^-3 (e.g., 1e17 cm^-3)
    n = 1e23
    # Electron mobility in m^2/(V*s) (e.g., 8000 cm^2/(V*s))
    mu = 0.8
    # Phase coherence time in seconds (strongly temperature-dependent, longer at low T)
    tau_phi = 1.5e-10 # (150 ps)

    # --- 3. Calculate Derived Physical Quantities ---
    # Effective mass
    m_star = m_star_ratio * m_e
    
    # Fermi velocity for a 3D degenerate electron gas
    v_F = (hbar / m_star) * (3 * math.pi**2 * n)**(1/3)
    
    # Momentum relaxation time (elastic scattering time) from mobility
    tau = (mu * m_star) / e
    
    # Elastic mean free path
    l = v_F * tau
    
    # Diffusion constant
    D = (1/3) * v_F**2 * tau
    
    # Phase coherence length
    L_phi = math.sqrt(D * tau_phi)

    # --- 4. Evaluate the Quantum Correction Formula ---
    # The formula for the 3D weak localization correction is:
    # Delta_sigma = -(e^2 / (pi^2 * hbar)) * (1/l - (1/L_phi) * atan(L_phi/l))
    
    prefactor = e**2 / (math.pi**2 * hbar)
    term_l = 1 / l
    term_L_phi = (1 / L_phi) * math.atan(L_phi / l)
    
    delta_sigma = -prefactor * (term_l - term_L_phi)

    # --- 5. Print the Results ---
    print("--- Evaluating Quantum Correction to Conductivity (Weak Localization in 3D) ---")
    print("\n[1] Assumed Semiconductor Parameters:")
    print(f"  - Effective Mass (m*): {m_star:.3e} kg")
    print(f"  - Carrier Density (n): {n:.1e} m^-3")
    print(f"  - Mobility (mu): {mu:.2f} m^2/Vs")
    print(f"  - Phase Coherence Time (tau_phi): {tau_phi:.1e} s")
    
    print("\n[2] Calculated Intermediate Values:")
    print(f"  - Elastic Mean Free Path (l): {l:.3e} m")
    print(f"  - Phase Coherence Length (L_phi): {L_phi:.3e} m")

    print("\n[3] Final Equation with Numerical Values:")
    print("  Δσ = - (e² / (π² * ħ)) * [1/l - (1/L_phi) * atan(L_phi/l)]")
    print(f"  Δσ = - ({e**2:.4e} / (π² * {hbar:.4e})) * [1/{l:.4e} - (1/{L_phi:.4e}) * atan({L_phi:.4e}/{l:.4e})]")
    print(f"  Δσ = - ({prefactor:.4e}) * [{term_l:.4e} - {term_L_phi:.4e}]")

    print("\n[4] Final Result:")
    print(f"The quantum correction to conductivity Δσ is: {delta_sigma:.4f} S/m (Siemens per meter).")

    # To fulfill the final output format requirement
    print(f"\n<<< {delta_sigma:.4f} >>>")

if __name__ == '__main__':
    evaluate_quantum_correction()
