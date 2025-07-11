import math

def calculate_conductivity_correction():
    """
    Calculates the quantum correction to conductivity (Weak Localization) for a
    bulk (3D) semiconductor.
    """
    # --- Physical Constants (SI units) ---
    e = 1.60217663e-19      # Elementary charge in Coulombs
    hbar = 1.05457182e-34   # Reduced Planck constant in J·s

    # --- Assumed System Parameters (typical for a doped semiconductor at low T) ---
    # l: elastic mean free path in meters.
    # Typical values are in the tens of nanometers.
    l = 50e-9  # 50 nm

    # L_phi: phase coherence length in meters.
    # Must be greater than l. It is strongly temperature-dependent.
    L_phi = 500e-9 # 500 nm

    print("This script evaluates the quantum correction to conductivity due to weak localization in a bulk semiconductor.")
    print("The correction is typically negative, representing a reduction in conductivity.")
    
    print("\n--- Input Parameters ---")
    print("Physical Constants:")
    print(f"  Elementary charge (e)       = {e:.4g} C")
    print(f"  Reduced Planck constant (ħ) = {hbar:.4g} J·s")
    print("\nAssumed Material Properties:")
    print(f"  Elastic mean free path (l)    = {l*1e9:.0f} nm = {l:.1e} m")
    print(f"  Phase coherence length (L_ϕ)  = {L_phi*1e9:.0f} nm = {L_phi:.1e} m")
    
    # --- Calculation ---
    # The formula is: Δσ = - (e² / (2π²ħ)) * (1/l - 1/L_ϕ)
    
    prefactor = -(e**2) / (2 * math.pi**2 * hbar)
    length_term = (1 / l) - (1 / L_phi)
    delta_sigma = prefactor * length_term

    # --- Output Results ---
    print("\n--- Calculation Steps ---")
    print("The formula is: Δσ = - (e² / (2π²ħ)) * (1/l - 1/L_ϕ)")
    
    print("\nFinal Equation with numerical values substituted:")
    print(f"Δσ = - (({e:.4g})² / (2 * ({math.pi:.4f})² * {hbar:.4g})) * (1 / {l:.1e} - 1 / {L_phi:.1e})")

    print("\n--- Final Result ---")
    print(f"The quantum correction to conductivity is: Δσ = {delta_sigma:.4f} S/m")

if __name__ == "__main__":
    calculate_conductivity_correction()
<<<Δσ = -222.0742 S/m>>>