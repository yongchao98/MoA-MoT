from fractions import Fraction

def solve_flux():
    """
    Calculates the energy flow through the yellow sides of the pyramid.
    """
    print("### Step-by-step Calculation of Energy Flow ###\n")

    # The problem is to find the flux of F = (3x^3*y^2*z, 3*x^2*y^3, z)
    # through the two yellow sides of a pyramid.
    # The pyramid has base vertices at (±1, ±1, 0) and an apex at (0, 0, 4).
    # "Interspersed" coloring means opposite faces have the same color.
    # We'll calculate the flux for the front (+y) and back (-y) faces.

    # --- Calculation for the Front Face (S1) ---
    # The plane is y = 1 - z/4. The normal dS is (0, 1, 1/4) dx dz.
    # The dot product F . dS = 3*x^2*y^3 + z/4.
    # Substituting y = 1 - z/4, the integrand is 3*x^2*(1 - z/4)^3 + z/4.
    # We integrate this over the projection on the xz-plane, which is a triangle
    # with vertices (-1, 0), (1, 0), and (0, 4).

    print("1. Calculating flux for the Front Face (S1):")

    # The integral splits into two parts.
    # Part A: integral of 3*x^2*(1 - z/4)^3
    # Analytically, this evaluates to 8/7.
    flux1_A = Fraction(8, 7)
    print(f"   - The result of the first part of the integral is {flux1_A.numerator}/{flux1_A.denominator}")

    # Part B: integral of z/4
    # Analytically, this evaluates to 4/3.
    flux1_B = Fraction(4, 3)
    print(f"   - The result of the second part of the integral is {flux1_B.numerator}/{flux1_B.denominator}")

    # Total flux for S1
    flux1 = flux1_A + flux1_B
    print(f"   - Total flux for the front face is the sum: {flux1_A.numerator}/{flux1_A.denominator} + {flux1_B.numerator}/{flux1_B.denominator} = {flux1.numerator}/{flux1.denominator}\n")

    # --- Calculation for the Back Face (S2) ---
    # The integrand and the domain of integration for the back face
    # are identical to the front face due to the symmetries of the
    # vector field and the surface geometry.
    flux2 = flux1
    print("2. Calculating flux for the Back Face (S2):")
    print(f"   - By symmetry, the flux through the back face is the same as the front face.")
    print(f"   - Total flux for the back face is {flux2.numerator}/{flux2.denominator}\n")


    # --- Total Flux ---
    total_flux = flux1 + flux2
    print("3. Calculating Total Flux through Yellow Sides:")
    print("   The total energy flow is the sum of the fluxes through the two yellow faces (front and back).")
    # Using the requirement: output each number in the final equation!
    print(f"\nFinal Equation: {flux1.numerator}/{flux1.denominator} + {flux2.numerator}/{flux2.denominator} = {total_flux.numerator}/{total_flux.denominator}")

    print(f"\nThe total energy flow is {total_flux.numerator}/{total_flux.denominator}, which is approximately {float(total_flux):.4f}.")

solve_flux()
<<<104/21>>>