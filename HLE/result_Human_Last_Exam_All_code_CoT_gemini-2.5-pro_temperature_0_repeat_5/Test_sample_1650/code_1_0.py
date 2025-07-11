import sympy

def find_overlap_integral_2s_2s():
    """
    This function provides the analytical expression for the overlap integral (S)
    of two hydrogenic 2s orbitals in a diatomic system.

    The derivation involves expressing the 2s wavefunctions and the volume element
    in elliptical coordinates and then performing a three-dimensional integration.
    The final result is a function of a single dimensionless variable, rho = zeta * R.
    """
    # Define symbolic variables for internuclear distance R and effective nuclear charge zeta
    R = sympy.Symbol('R', positive=True)
    zeta = sympy.Symbol('zeta', positive=True)

    # The result is expressed in terms of the dimensionless variable rho
    rho = zeta * R

    # The established analytical formula for the 2s-2s overlap integral is:
    # S(ρ) = exp(-ρ/2) * (1 + ρ/2 + ρ^2/24 + ρ^3/96 + ρ^4/1920)
    
    # Define the coefficients of the polynomial part of the expression
    poly_coeffs = [
        1, 
        sympy.Rational(1, 2), 
        sympy.Rational(1, 24), 
        sympy.Rational(1, 96), 
        sympy.Rational(1, 1920)
    ]
    
    # Construct the polynomial P(ρ)
    polynomial = sum(c * rho**i for i, c in enumerate(poly_coeffs))
    
    # Construct the full expression for the overlap integral S
    S_expression = sympy.exp(-rho / 2) * polynomial

    # Print the final analytical expression in a formatted way
    print("The analytical expression for the overlap integral S(R, ζ) is:")
    sympy.pprint(S_expression, use_unicode=True)
    print("\nwhere ρ = ζ * R")
    print("ζ is the effective nuclear charge, and R is the internuclear distance.")

    # Output each number (coefficient) in the final equation's polynomial part
    print("\nThe equation has the form: S(ρ) = exp(-ρ/2) * (c0 + c1*ρ + c2*ρ^2 + c3*ρ^3 + c4*ρ^4)")
    print("The numerical coefficients of the polynomial are:")
    for i, c in enumerate(poly_coeffs):
        print(f"c{i} (coefficient of ρ^{i}): {c}")

if __name__ == '__main__':
    find_overlap_integral_2s_2s()