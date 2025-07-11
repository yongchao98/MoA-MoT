def solve_edl_potential():
    """
    This function derives and prints the expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y) in a parallel-plate microchannel under specific
    boundary conditions.
    """

    # --- Introduction and Problem Setup ---
    print("### Derivation of the Electrical Double-Layer (EDL) Potential Distribution psi(y) ###")
    print("-" * 75)
    print("This script solves for the EDL potential psi(y) based on the following:")
    print("\n1. Governing Equation (Linearized Poisson-Boltzmann Equation):")
    print("   d^2(psi)/dy^2 = k^2 * psi")
    print("\n   where:")
    print("   - psi: Electrical potential")
    print("   - y:   Coordinate perpendicular to the plates (from -H/2 to H/2)")
    print("   - k:   Debye-Huckel parameter")
    print("   - H:   Channel height")

    print("\n2. Boundary Conditions:")
    print("   The potential at the plates is given by the slip-dependent zeta potential.")
    print("   - Top plate (y = H/2):")
    print("     psi(H/2) = z_a2 = z_2 * (1 + beta * k)")
    print("     Given z_2 = 0, the boundary condition becomes: psi(H/2) = 0")
    print("\n   - Bottom plate (y = -H/2):")
    print("     psi(-H/2) = z_a1 = z_1 * (1 + beta * k)")

    print("\n3. Solving the Differential Equation:")
    print("   The general solution to the governing equation is:")
    print("   psi(y) = C*cosh(k*y) + D*sinh(k*y)")
    print("\n   Applying the two boundary conditions allows us to solve for the constants C and D.")
    print("   The result is a specific solution for psi(y).")
    print("-" * 75)

    # --- Final Expression ---
    print("\n### Final Expression for EDL Potential Distribution ###")
    print("-" * 75)
    
    # Define the components of the final equation as strings
    z_a1 = "z_1 * (1 + beta * k)"
    numerator = "sinh(k * (H/2 - y))"
    denominator = "sinh(k * H)"
    
    # Print the final expression step-by-step
    print(f"The potential distribution psi(y) is given by:")
    print(f"\n   psi(y) = [ {z_a1} ] * [ {numerator} / {denominator} ]\n")

    # --- Breakdown of the Final Equation ---
    print("Breakdown of the expression's components:")
    print("------------------------------------------")
    
    term1_desc = "The effective zeta potential at the bottom wall (y = -H/2):"
    term1_val = "z_a1 = z_1 * (1 + beta * k)"
    print(f"- {term1_desc:<60} {term1_val}")
    
    term2_desc = "Numerator term, describing the decay from the top wall:"
    term2_val = numerator
    print(f"- {term2_desc:<60} {term2_val}")

    term3_desc = "Denominator term, for normalization:"
    term3_val = denominator
    print(f"- {term3_desc:<60} {term3_val}")

if __name__ == '__main__':
    solve_edl_potential()
