def solve_mutual_inductance_change():
    """
    This function calculates and prints the expression for the change in mutual
    inductance per unit length (ΔM') between two wire pairs when surrounded
    by magnetic concentrators.
    """

    # Explain the context of the calculation
    print("The expression for the change in mutual inductance per unit length (ΔM') is derived using the dipole approximation (d >> h).")
    print("ΔM' is the difference between the inductance with concentrators (M₂') and the bare inductance (M₁').")
    print("-" * 50)

    # Define the components of the formula using Unicode for clarity
    delta_M_prime = "ΔM'"
    h_sq = "h²"
    d_sq = "d²"
    R2_div_R1 = "R₂/R₁"

    # The simplified coefficient C = μ₀ / (2π) = (4π × 10⁻⁷) / (2π) = 2 × 10⁻⁷ H/m
    coefficient = "2 × 10⁻⁷"

    # Assemble and print the final equation.
    # This format explicitly shows each "number" and variable in the equation.
    print("The final expression for the change in mutual inductance per unit length is:")
    print(f"{delta_M_prime} = ( {coefficient} * {h_sq} / {d_sq} ) * ( {R2_div_R1} - 1 )")

    print("\nWhere:")
    print("  ΔM' = The change in mutual inductance per unit length (in Henries/meter)")
    print(f"  {coefficient} = The magnetic constant μ₀/2π (in Henries/meter)")
    print("  h = The separation between the wires in each circuit (in meters)")
    print("  d = The distance between the centers of the two circuits (in meters)")
    print("  R₁ = The inner radius of the concentrator shells (in meters)")
    print("  R₂ = The outer radius of the concentrator shells (in meters)")

solve_mutual_inductance_change()