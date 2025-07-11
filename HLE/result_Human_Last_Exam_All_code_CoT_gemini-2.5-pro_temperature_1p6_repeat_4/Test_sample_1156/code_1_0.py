import math

def solve():
    """
    This function calculates the normalization constant for the invariant density
    and prints the final equation.
    """
    # The analytical form of the density is rho(x) = C / (1 + 2*sqrt(x)).
    # The normalization constant C is found by integrating rho(x) from 0 to 1
    # and setting the result to 1.
    # C * integral(1/(1+2*sqrt(x))) dx from 0 to 1 = 1
    # The value of the integral is 1 - ln(3)/2.
    # So, C = 1 / (1 - ln(3)/2) = 2 / (2 - ln(3)).

    # Constants in the equation for rho(x)
    numerator_C = 2
    denominator_C_term1 = 2
    denominator_C_term2 = 3 # Inside the ln
    denominator_rho_term1 = 1
    denominator_rho_term2 = 2

    print("The normalised density of the invariant measure is:")
    print(f"rho(x) = ({numerator_C} / ({denominator_C_term1} - ln({denominator_C_term2}))) / ({denominator_rho_term1} + {denominator_rho_term2}*sqrt(x))")

solve()