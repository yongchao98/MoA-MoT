import math

def solve_shock_profile():
    """
    Calculates and prints the analytical solution for the density profile of a
    Mach sqrt(3) shock wave in a monatomic ideal gas with Pr = 3/4.
    """

    # Define constants from the problem statement
    # Mach number, M
    M = math.sqrt(3)
    # Adiabatic index for a monatomic gas, gamma
    gamma = 5.0 / 3.0

    # The analytical solution for the normalized density (rho/rho_0) as a function
    # of normalized position (x/L) is found to be of the form:
    # rho/rho_0 = 1 + sqrt(exp(coeff * c_0 * x/L) / (3 + exp(coeff * c_0 * x/L)))
    # where c_0 is the ambient sound speed.

    # The coefficient in the exponent is derived from the governing equations.
    # Its value is (M * (gamma + 1)) / 4.
    coeff = M * (gamma + 1) / 4.0

    # As requested, output each number in the final equation.
    print("The numbers that form the final analytical equation are:")
    print(f"The additive constant: 1")
    print(f"The constant in the denominator of the fraction under the square root: 3")
    print(f"The numerical coefficient in the exponent: {coeff}")
    print("-" * 50)

    # Print the final equation for the density profile.
    # rho_norm represents rho/rho_0
    # x_norm represents x/L
    print("The analytical solution for the density profile is:")
    print("rho/rho_0 = 1 + (exp(coeff * c_0 * x/L) / (3 + exp(coeff * c_0 * x/L)))**0.5")
    print("\nSubstituting the value of the coefficient (coeff):")
    print(f"rho/rho_0 = 1 + (exp({coeff:.5f} * c_0 * x/L) / (3 + exp({coeff:.5f} * c_0 * x/L)))**0.5")
    print("\nWhere:")
    print("  - rho/rho_0 is the density normalized by the ambient density.")
    print("  - x/L is the position normalized by the ambient conductive length scale.")
    print("  - c_0 is the ambient sound speed.")
    print("  - exp() is the exponential function, and **0.5 denotes the square root.")

solve_shock_profile()