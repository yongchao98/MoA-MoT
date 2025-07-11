import numpy as np
from scipy.constants import mu_0

def calculate_inductance_change():
    """
    Calculates the change in mutual inductance between two circuits when
    magnetic concentrators are added.
    """
    # --- Define System Parameters ---
    # h: separation of wires in each circuit (m)
    h = 0.01
    # d: distance between the two circuits (m)
    d = 0.5
    # R1: inner radius of the concentrator shell (m)
    R1 = 0.02
    # R2: outer radius of the concentrator shell (m)
    R2 = 0.04
    # Note: The calculation is performed per unit length, as the wires are "very long".

    # --- Derivation Summary ---
    # 1. The mutual inductance per unit length for the bare circuits (M₁') in the
    #    far-field approximation (d >> h) is: M₁' = (μ₀ * h²) / (2 * π * d²)
    # 2. Each ideal concentrator enhances the magnetic coupling by a factor of (R₂/R₁).
    #    With two concentrators (one on each circuit), the total enhancement to the
    #    mutual inductance is (R₂/R₁)²
    # 3. The new mutual inductance (M₂') is: M₂' = M₁' * (R₂/R₁)²
    # 4. The change in mutual inductance (ΔM') is M₂' - M₁', which simplifies to:
    #    ΔM' = M₁' * [(R₂/R₁)² - 1]

    # --- Calculation ---
    try:
        # Calculate the bare mutual inductance per unit length
        m1_per_L = (mu_0 * h**2) / (2 * np.pi * d**2)

        # Calculate the ratio of the radii
        radii_ratio_sq = (R2 / R1)**2

        # Calculate the change in mutual inductance per unit length
        delta_M_per_L = m1_per_L * (radii_ratio_sq - 1)

        # --- Output the results ---
        print("The expression for the change in mutual inductance per unit length (ΔM') is:")
        print("ΔM' = (μ₀ * h² / (2 * π * d²)) * [ (R₂/R₁)² - 1 ]")
        print("\nSubstituting the given values into the equation:")
        
        # Constructing the equation string with numbers
        # This string shows each number used in the final calculation
        equation_str = (f"ΔM' = ({mu_0:.3e} H/m * ({h} m)² / (2 * {np.pi:.4f} * ({d} m)²)) "
                        f"* [ ({R2} m / {R1} m)² - 1 ]")
        
        print(equation_str)
        print(f"\nCalculated change in mutual inductance per unit length:")
        print(f"ΔM' = {delta_M_per_L:.3e} H/m")

    except ZeroDivisionError:
        print("Error: Distances and radii (d, R1) cannot be zero.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Execute the calculation
calculate_inductance_change()