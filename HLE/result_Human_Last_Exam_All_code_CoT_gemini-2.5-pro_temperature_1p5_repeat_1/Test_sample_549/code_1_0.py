import numpy as np
import scipy.constants as const

def evaluate_quantum_correction():
    """
    Calculates the quantum correction to conductivity in a 3D semiconductor
    based on the weak localization theory.
    """
    # --- Step 1: Define constants and parameters ---
    # Physical constants from scipy.constants
    e = const.e  # Elementary charge in Coulombs
    hbar = const.hbar  # Reduced Planck constant in J.s
    pi = np.pi

    # Typical parameters for a doped semiconductor at low temperature.
    # These values can vary significantly based on material, doping, and temperature.
    l_e = 50e-9  # Elastic mean free path in meters (50 nm)
    L_phi = 500e-9 # Phase coherence length in meters (500 nm)

    # --- Step 2: Perform the calculation ---
    # The prefactor e^2 / (2 * π^2 * ħ) has units of Siemens (S) or (Ω)^-1
    prefactor = e**2 / (2 * pi**2 * hbar)
    
    # The term (1/L_φ - 1/lₑ) has units of m^-1
    length_term = (1 / L_phi) - (1 / l_e)
    
    # The final quantum correction Δσ has units of S/m
    delta_sigma = prefactor * length_term

    # --- Step 3: Print the results ---
    print("Evaluation of the Quantum Correction to Conductivity (Δσ)")
    print("---------------------------------------------------------")
    print("The formula for weak localization in 3D is:")
    print("Δσ = (e² / (2 * π² * ħ)) * (1/L_φ - 1/lₑ)\n")

    print("Using the following physical constants and typical parameters:")
    print(f"  e (elementary charge)   = {e:.5e} C")
    print(f"  ħ (reduced Planck const)  = {hbar:.5e} J·s")
    print(f"  lₑ (elastic mean free path) = {l_e:.1e} m ({l_e*1e9:.0f} nm)")
    print(f"  L_φ (phase coherence length) = {L_phi:.1e} m ({L_phi*1e9:.0f} nm)\n")

    # Output the equation with numbers plugged in, as requested.
    print("Substituting the values into the formula:")
    print(f"Δσ = ({e**2:.5e} C² / (2 * {pi**2:.4f} * {hbar:.5e} J·s)) * (1/{L_phi:.1e} m⁻¹ - 1/{l_e:.1e} m⁻¹)")
    print(f"Δσ = ({prefactor:.5e} S) * ({length_term:.3e} m⁻¹)")
    
    print("\n--- Result ---")
    print(f"The calculated quantum correction to conductivity is: Δσ = {delta_sigma:.2f} S/m")
    print("\nThe negative sign indicates a reduction in conductivity, which is characteristic of the weak localization effect.")

if __name__ == '__main__':
    evaluate_quantum_correction()

<<< -222.00 >>>