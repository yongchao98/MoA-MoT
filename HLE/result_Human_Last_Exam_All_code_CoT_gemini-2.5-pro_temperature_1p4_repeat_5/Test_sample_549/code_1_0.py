import numpy as np

def calculate_conductivity_correction():
    """
    Calculates the quantum correction to conductivity for an electron in a 3D
    disordered medium (weak localization).
    """

    # Physical constants in SI units
    e = 1.60217663e-19  # Elementary charge in Coulombs
    hbar = 1.05457182e-34 # Reduced Planck constant in J·s

    # Typical parameters for a bulk semiconductor at low temperature
    # Elastic mean free path (l) in meters
    l = 100e-9
    # Phase coherence length (L_phi) in meters
    L_phi = 1000e-9

    # The formula for the 3D quantum correction to conductivity is:
    # Δσ = [e² / (2π²ħ)] * (1/L_phi - 1/l)

    # Calculate the prefactor
    prefactor = e**2 / (2 * np.pi**2 * hbar)

    # Calculate the length-dependent term
    length_term = (1/L_phi) - (1/l)

    # Calculate the final quantum correction
    delta_sigma = prefactor * length_term

    # --- Outputting the step-by-step calculation ---
    print("Evaluating the quantum correction to conductivity (Δσ) in a bulk semiconductor.")
    print("\n1. The formula for the 3D weak localization correction is:")
    print("   Δσ = [e² / (2π²ħ)] * (1/L_phi - 1/l)")

    print("\n2. Plugging in the values:")
    print(f"   e (elementary charge) = {e:.4e} C")
    print(f"   ħ (reduced Planck constant) = {hbar:.4e} J·s")
    print(f"   l (mean free path) = {l:.0e} m")
    print(f"   L_phi (phase coherence length) = {L_phi:.0e} m")

    print("\n3. Calculating each part of the equation:")
    print(f"   Prefactor [e² / (2π²ħ)] = ({e:.4e})² / (2 * π² * {hbar:.4e}) = {prefactor:.4e} S")
    print(f"   Length term (1/L_phi - 1/l) = (1/{L_phi:.0e} - 1/{l:.0e}) = {length_term:.2e} m⁻¹")
    
    print("\n4. Final Calculation:")
    print(f"   Δσ = {prefactor:.4e} S * ({length_term:.2e} m⁻¹)")
    print(f"   Δσ = {delta_sigma:.2f} S/m")
    
    return delta_sigma

# Run the calculation and store the result
final_answer = calculate_conductivity_correction()
# The final answer needs to be enclosed in <<<>>>
# print(f"<<<{final_answer:.1f}>>>")