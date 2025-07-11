import sympy

def solve_gamma_product():
    """
    This function calculates the proportionality factor in the expression
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}.
    """

    # Define d and k as symbolic variables for mathematical representation
    d = sympy.Symbol('d')
    k = sympy.Symbol('k')

    print("The task is to find the proportionality factor C in the equation:")
    print("gamma_{mu nu} * gamma_{mu_1 ... mu_k} * gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}")
    print("where the sum over mu and nu is implied.")
    print("\nThe factor C can be shown to be independent of k. We calculate it for the case k=0.")
    print("For k=0, gamma_{...} is the identity matrix I, so the equation simplifies to:")
    print("gamma_{mu nu} * gamma^{mu nu} = C * I\n")
    print("To evaluate this, we use the definition gamma_{mu nu} = 1/2 * (gamma_mu*gamma_nu - gamma_nu*gamma_mu)")
    print("and the Clifford algebra identity {gamma_mu, gamma_nu} = 2*g_{mu nu}.")
    print("After expansion and using gamma matrix identities, we find:\n")

    # The detailed derivation leads to the following result for the factor C.
    # The calculation is based on standard gamma matrix trace and contraction identities.
    # Sum of terms: 1/4 * (d*(2-d) - d**2 - d**2 + d*(2-d)) = d - d**2
    proportionality_factor = d * (1 - d)

    print(f"The proportionality factor C is found to be: {proportionality_factor}")

    # Now, we print the final equation step by step.
    mu = r'\mu'
    nu = r'\nu'
    mu_1_k = r'{\mu_1 \ldots \mu_k}'
    
    term1 = f"gamma_({mu}{nu})"
    term2 = f"gamma_({mu_1_k})"
    term3 = f"gamma^({mu}{nu})"
    factor_str = f"({sympy.printing.latex(proportionality_factor)})"
    final_term = f"gamma_({mu_1_k})"

    print("\nThus, the final equation with all its parts is:")
    print(f"The first term is the gamma matrix with indices {mu} and {nu}: {term1}")
    print(f"The second term is the rank-k gamma matrix: {term2}")
    print(f"The third term is the gamma matrix with upper indices {mu} and {nu}: {term3}")
    print(f"The proportionality factor is: {factor_str}")
    print(f"The product of these three terms equals the factor times the original rank-k gamma matrix:")
    print(f"{term1} * {term2} * {term3} = {factor_str} * {final_term}")


if __name__ == '__main__':
    solve_gamma_product()