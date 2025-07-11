import math

def solve_mutual_inductance_change():
    """
    This script derives and prints the expression for the change in mutual inductance (ΔM)
    between two circuits when placed inside magnetic concentrators.
    The result is per unit length and valid in the limit where d >> h.
    """

    # --- Symbolic Representation of Physical Quantities ---
    # These strings represent the variables in our final equation.
    mu_0 = "μ₀"  # Permeability of free space
    pi = "π"      # Mathematical constant Pi
    h = "h"       # Wire separation in each circuit
    d = "d"       # Distance between the two circuits
    R1 = "R₁"     # Inner radius of the concentrator shells
    R2 = "R₂"     # Outer radius of the concentrator shells

    print("Step 1: The mutual inductance per unit length for the bare circuits (M₁) is derived.")
    # In the limit d >> h, M₁/L ≈ (μ₀ * h²) / (2 * π * d²)
    m1_factor = f"({mu_0} * {h}**2) / (2 * {pi} * {d}**2)"
    print(f"M₁/L ≈ {m1_factor}\n")

    print("Step 2: The effect of the concentrators is to scale the effective wire separation.")
    # The scaling factor is (R₂/R₁)
    enhancement_factor = f"({R2}/{R1})"
    print(f"Effective wire separation h_eff = {h} * {enhancement_factor}\n")

    print("Step 3: The mutual inductance with concentrators (M₂) is calculated.")
    # M₂ is proportional to h_eff², so it's scaled by (R₂/R₁)²
    m2_expression = f"{m1_factor} * {enhancement_factor}**2"
    print(f"M₂/L ≈ {m2_expression}\n")

    print("Step 4: The change in mutual inductance (ΔM = M₂ - M₁) is calculated.")
    # ΔM/L = M₂/L - M₁/L = M₁/L * [ (R₂/R₁)² - 1 ]
    print("The final expression for the change in mutual inductance per unit length (ΔM/L) is:\n")
    
    # --- Final Equation Output ---
    # The instruction "output each number in the final equation" is interpreted as
    # showing the individual components and numbers that make up the formula.
    
    print("Component parts of the equation:")
    print(f"  - Bare Inductance Term: {m1_factor}")
    print(f"  - Concentrator Factor: (({R2}/{R1})**2 - 1)")
    
    print("\nLiteral numbers in the equation:")
    print("  - Number: 2 (in the denominator of the first term)")
    print("  - Number: 2 (in the exponent of the concentrator factor)")
    print("  - Number: 1 (subtracted in the concentrator factor)")
    
    print("\n-------------------------------------------------------------")
    print("Final Assembled Equation for ΔM (per unit length):")
    # This prints the final complete expression as the script's main output.
    final_equation = f"{m1_factor} * (({R2}/{R1})**2 - 1)"
    print(final_equation)
    print("-------------------------------------------------------------")

solve_mutual_inductance_change()