import math

def calculate_conductivity_correction():
    """
    Calculates and prints the quantum correction to conductivity for an electron
    in a bulk semiconductor based on the theory of weak localization.
    """
    # Physical Constants in SI units
    e = 1.602176634e-19      # Elementary charge in Coulombs (C)
    hbar = 1.054571817e-34   # Reduced Planck constant in Joule-seconds (J·s)
    pi = math.pi

    # Typical parameters for a doped semiconductor at low temperature
    # l: elastic mean free path in meters (m)
    l = 50e-9
    # L_phi: phase coherence length in meters (m)
    L_phi = 500e-9

    # --- Calculation ---
    # The prefactor e^2 / (2 * pi^2 * hbar) has units of Siemens (S)
    prefactor = e**2 / (2 * pi**2 * hbar)
    # The term (1/l - 1/L_phi) has units of inverse meters (m⁻¹)
    length_term = (1/l - 1/L_phi)
    # The final correction to conductivity (delta_sigma) is in S/m
    delta_sigma = -prefactor * length_term

    # --- Output ---
    print("This script evaluates the quantum correction to conductivity (Δσ) in a 3D bulk semiconductor.")
    print("The correction is due to weak localization, a quantum interference effect.\n")
    print("The governing equation is:")
    print("  Δσ = - (e² / (2 * π² * ħ)) * (1/l - 1/L_φ)\n")

    print("Using the following values:")
    print(f"  e (Charge)             = {e:.4e} C")
    print(f"  ħ (Reduced Planck)     = {hbar:.4e} J·s")
    print(f"  l (Mean Free Path)     = {l:.1e} m")
    print(f"  L_φ (Phase Coherence)  = {L_phi:.1e} m\n")

    print("Substituting the numbers into the equation:")
    print(f"  Δσ = - (({e:.4e})² / (2 * π² * {hbar:.4e})) * (1/{l:.1e} - 1/{L_phi:.1e})")
    print(f"  Δσ = - ({prefactor:.4e} S) * ({length_term:.4e} m⁻¹)")
    
    print("\n---------------------------------------------------")
    print(f"Calculated Quantum Correction: Δσ = {delta_sigma:.2f} S/m")
    print("---------------------------------------------------")
    print("\nThis negative value indicates a reduction in conductivity (increase in resistivity) due to quantum effects.")

if __name__ == '__main__':
    calculate_conductivity_correction()
    # The final answer format is handled outside the function for clarity.
    # Calculation result: delta_sigma = -1.2339e-5 * 1.8e7 = -222.10
    final_answer = -222.10

    # The prompt asks for the final answer in a specific format.
    # To conform, we will output it here as requested.
    # print(f"\n<<<{final_answer:.1f}>>>") # Example for display