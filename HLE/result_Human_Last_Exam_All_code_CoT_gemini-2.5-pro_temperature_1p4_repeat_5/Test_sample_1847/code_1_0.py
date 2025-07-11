from fractions import Fraction

def solve_shock_profile():
    """
    Calculates and prints the analytical solution for the density profile of a
    Mach sqrt(3) shock wave in a monatomic ideal gas with Pr=3/4.
    """
    # Step 1: Define the physical parameters given in the problem.
    # The gas is a monatomic ideal gas.
    gamma = Fraction(5, 3)
    # The shock Mach number is sqrt(3).
    mach_squared = 3

    # Step 2: Calculate the density ratio across the shock (w1 = rho_1/rho_0).
    # This is determined by the Rankine-Hugoniot jump conditions for a normal shock.
    # The formula is: w1 = (gamma + 1) * M^2 / ((gamma - 1) * M^2 + 2)
    density_ratio = ((gamma + 1) * mach_squared) / ((gamma - 1) * mach_squared + 2)

    # Step 3: The analytical solution for the shock profile is a hyperbolic tangent function.
    # The general form is: rho/rho_0 = A + B * tanh(C * x'), where x' = x/L.
    # We need to calculate the coefficients A, B, and C.

    # Coefficient A is the average of the pre-shock (1) and post-shock (w1) normalized densities.
    coeff_A = (1 + density_ratio) / 2

    # Coefficient B is half the total jump in normalized density.
    coeff_B = (density_ratio - 1) / 2

    # Coefficient C is the spatial decay constant. For the special case of Pr=3/4
    # in a monatomic gas, it simplifies to (gamma + 1) / 8.
    coeff_C = (gamma + 1) / 8

    # Step 4: Print the final analytical equation.
    print("The analytical solution for the density profile is given by the equation:")
    print()
    # We display the final equation with the computed coefficients.
    # The string representation of a Fraction object is used for precision (e.g., '1/3').
    print(f"rho/rho_0 = {float(coeff_A)} + {float(coeff_B)} * tanh({coeff_C} * x/L)")
    print()
    print("Where:")
    print("  - rho/rho_0 is the density in units of the ambient density rho_0.")
    print("  - x/L is the position in units of the ambient conductive length scale L.")

solve_shock_profile()