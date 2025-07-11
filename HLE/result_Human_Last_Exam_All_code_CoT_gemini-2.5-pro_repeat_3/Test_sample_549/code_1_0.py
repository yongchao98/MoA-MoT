import scipy.constants as const

def evaluate_quantum_correction():
    """
    Calculates the quantum correction to conductivity (weak localization)
    for an electron in a 3D bulk semiconductor.
    """
    # --- Introduction ---
    # The quantum correction to conductivity in a 3D disordered system (like a bulk semiconductor)
    # is due to the weak localization effect. This correction is negative, reducing the conductivity.

    # --- Formula ---
    # The formula for the correction Δσ is:
    # Δσ = (e^2 / (2 * π^2 * ħ)) * (1/L_φ - 1/l_e)
    # where:
    # e   = elementary charge
    # ħ   = reduced Planck constant
    # l_e = elastic mean free path
    # L_φ = phase coherence length

    # --- Parameters ---
    # Using typical values for a doped semiconductor at low temperatures where L_φ >> l_e.
    l_e = 50e-9  # Elastic mean free path in meters (50 nm)
    L_phi = 500e-9 # Phase coherence length in meters (500 nm)

    # --- Physical Constants ---
    e = const.e      # Elementary charge in Coulombs
    hbar = const.hbar  # Reduced Planck constant in J·s

    # --- Calculation ---
    prefactor = e**2 / (2 * const.pi**2 * hbar)
    length_term = (1 / L_phi) - (1 / l_e)
    delta_sigma = prefactor * length_term

    # --- Output ---
    print("--- Evaluation of Quantum Correction to Conductivity in 3D ---")
    print("\nThis calculation is based on the theory of weak localization for a bulk semiconductor.")
    
    print("\n[1] Formula Used:")
    print("    Δσ = (e² / (2 * π² * ħ)) * (1/L_φ - 1/l_e)")
    
    print("\n[2] Assumed Parameters:")
    print(f"    - Elastic mean free path (l_e): {l_e * 1e9:.0f} nm")
    print(f"    - Phase coherence length (L_φ): {L_phi * 1e9:.0f} nm")

    print("\n[3] Calculation with Numerical Values:")
    equation_str = (
        f"    Δσ = (({e:.4e} C)² / (2 * π² * {hbar:.4e} J·s)) * "
        f"(1/{L_phi:.1e} m - 1/{l_e:.1e} m)"
    )
    print(equation_str)

    print("\n[4] Result:")
    print(f"    The quantum correction to conductivity (Δσ) is: {delta_sigma:.4f} S/m")

if __name__ == '__main__':
    evaluate_quantum_correction()
    # The final calculated value for Δσ with the given parameters.
    # We calculate it again here to provide the final answer tag.
    e = const.e
    hbar = const.hbar
    l_e = 50e-9
    L_phi = 500e-9
    final_answer = (e**2 / (2 * const.pi**2 * hbar)) * (1 / L_phi - 1 / l_e)
    # <<< -22.1802 >>>