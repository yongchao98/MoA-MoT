import math

def solve_shock_profile():
    """
    This function provides the analytical solution for the density profile of a
    Mach sqrt(3) shock wave in a monatomic ideal gas with Pr = 3/4.
    The code constructs and prints the final equation in a formatted string.
    """

    # The analytical derivation for the given conditions leads to the following
    # relationship between the normalized density (ρ/ρ₀) and the normalized
    # position (x/L).

    # The numbers that appear in the final equation are defined here.
    # The equation form is: ρ/ρ₀ = c1 + c2 / sqrt(c3 + c4 * exp( (c5_num/c5_den) * x/L ) )
    c1 = 1
    c2 = 1
    c3 = 1
    c4 = 3
    c5_num = -4
    c5_den = 3

    # Construct the formatted string for the final equation and explanation.
    equation_string = (
        "The analytical solution for the normalized density profile ρ/ρ₀ as a function of\n"
        "the normalized position x/L is:\n\n"
        f"ρ/ρ₀ = {c1} + {c2} / sqrt( {c3} + {c4} * exp( ({c5_num}/{c5_den}) * x/L ) )\n\n"
        "Where:\n"
        "  ρ   is the local density.\n"
        "  ρ₀  is the ambient (pre-shock) density.\n"
        "  x   is the position along the direction of propagation.\n"
        "  L   is the ambient conductive length scale, L = κ/(ρ₀*M*c₀*Cᵥ).\n"
        "  sqrt() is the square root function.\n"
        "  exp() is the exponential function.\n\n"
        "The numbers in this final equation are:\n"
        f"  Additive constant: {c1}\n"
        f"  Numerator in the fraction: {c2}\n"
        f"  Constant term inside the square root: {c3}\n"
        f"  Coefficient of the exponential term: {c4}\n"
        f"  Coefficient in the exponent: {c5_num}/{c5_den}"
    )

    print(equation_string)

solve_shock_profile()
<<<ρ/ρ₀ = 1 + 1 / sqrt( 1 + 3 * exp( (-4/3) * x/L ) )>>>