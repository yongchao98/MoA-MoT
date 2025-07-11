import math

def calculate_2s_2s_overlap():
    """
    Calculates the overlap integral for two 2s orbitals in H2+.
    
    This function uses the analytical formula derived using elliptical coordinates
    and prints the step-by-step calculation.
    """
    # Parameters for the H2+ ion in atomic units
    # R: internuclear distance in Bohr radii (a0)
    # zeta: effective nuclear charge
    R = 1.997  # Equilibrium bond length of H2+ is ~1.997 a0
    zeta = 1.24  # A common variational parameter for H2+

    print(f"Calculating the 2s-2s overlap integral for R = {R} a0 and zeta = {zeta}\n")

    # Calculate the dimensionless parameter p
    p = (zeta * R) / 2.0

    print(f"The formula is S = exp(-p) * (1 + p + p^2/3 + p^4/15), where p = zeta * R / 2.")
    print(f"For the given values, p = ({zeta} * {R}) / 2 = {p:.4f}\n")

    # Calculate each term in the polynomial part of the formula
    term1 = 1.0
    term2 = p
    term3 = p**2 / 3.0
    term4 = p**4 / 15.0

    # Calculate the sum of the polynomial terms
    poly_sum = term1 + term2 + term3 + term4

    # Calculate the exponential part
    exp_part = math.exp(-p)

    # Calculate the final overlap integral
    overlap_integral = exp_part * poly_sum

    # Print the equation with the calculated numbers
    print("Plugging the numbers into the formula:")
    equation_str = (
        f"S = exp(-{p:.4f}) * (1 + {p:.4f} + {p:.4f}^2/3 + {p:.4f}^4/15)"
    )
    print(equation_str)
    
    # Print the equation with each term evaluated
    evaluated_terms_str = (
        f"S = {exp_part:.4f} * ({term1:.4f} + {term2:.4f} + {term3:.4f} + {term4:.4f})"
    )
    print(evaluated_terms_str)
    
    sum_str = f"S = {exp_part:.4f} * ({poly_sum:.4f})"
    print(sum_str)

    # Print the final result
    print(f"\nThe final overlap integral S is: {overlap_integral:.4f}")

# Execute the function
calculate_2s_2s_overlap()
<<<0.8424>>>