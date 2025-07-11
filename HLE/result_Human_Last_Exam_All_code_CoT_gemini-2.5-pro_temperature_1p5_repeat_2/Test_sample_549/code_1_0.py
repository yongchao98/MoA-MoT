import numpy as np

def evaluate_quantum_correction():
    """
    Calculates the weak localization correction to conductivity in a 3D system.
    """
    # Physical constants in SI units
    e = 1.60217663e-19      # Elementary charge in Coulombs
    hbar = 1.05457182e-34   # Reduced Planck constant in J·s

    # Illustrative parameters for a bulk semiconductor at low temperature
    # l: mean free path (distance between elastic scattering events)
    l = 50e-9  # meters (50 nm)
    # L_phi: phase coherence length (distance before phase randomization)
    L_phi = 1000e-9 # meters (1 µm)

    # The formula for the quantum correction to conductivity (Δσ) in 3D is:
    # Δσ = - [e² / (2 * π² * ħ)] * (1/l - 1/L_φ)

    # --- Calculation ---
    prefactor = e**2 / (2 * np.pi**2 * hbar)
    length_term = (1/l - 1/L_phi)
    delta_sigma = -prefactor * length_term

    # --- Output Results ---
    print("Evaluation of Quantum Correction to Conductivity (Weak Localization) in 3D\n")
    print("Formula: Δσ = - [e² / (2 * π² * ħ)] * (1/l - 1/L_φ)\n")
    
    print("Using the following physical constants and illustrative parameters:")
    print(f"  e (elementary charge)      = {e:.4g} C")
    print(f"  ħ (reduced Planck constant)  = {hbar:.4g} J·s")
    print(f"  l (mean free path)         = {l:.1e} m")
    print(f"  L_φ (phase coherence length) = {L_phi:.1e} m\n")

    print("Step-by-step substitution into the formula:")
    print(f"Δσ = - [({e:.4g} C)² / (2 * {np.pi:.5f}² * {hbar:.4g} J·s)] * [1 / {l:.1e} m - 1 / {L_phi:.1e} m]")
    
    print("\nIntermediate values:")
    print(f"   The prefactor [e² / (2 * π² * ħ)] = {prefactor:.4g} S (Siemens)")
    print(f"   The length term [(1/l - 1/L_φ)]  = {length_term:.4g} 1/m")

    print("\n----------------------------------------------------")
    print(f"Final Result for the Quantum Correction:")
    print(f"Δσ = {delta_sigma:.3f} S/m")
    print("----------------------------------------------------\n")
    print("This negative correction indicates that quantum interference slightly increases the material's resistivity.")

if __name__ == '__main__':
    evaluate_quantum_correction()
