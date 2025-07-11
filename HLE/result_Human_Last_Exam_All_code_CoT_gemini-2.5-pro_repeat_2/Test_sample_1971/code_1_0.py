import sympy

def solve_supergravity_parameters():
    """
    This function calculates and prints the values of beta and alpha^2
    based on the principles of supersymmetry in N=1, d=4 supergravity.
    """

    # 1. Determine beta
    # The cancellation of S-linear terms in the variation of the
    # super-cosmological constant leads to the equation:
    # 1/4 + beta = 0
    beta = -1/4

    # 2. Determine alpha^2
    # By analyzing the bosonic vacuum of the theory, we find the relation
    # between the Ricci scalar R and alpha:
    # R = -3 * kappa^2 * alpha^2
    # So, alpha^2 = -R / (3 * kappa^2)

    # The problem asks for the number of beta and the number of alpha^2 in terms of R.
    # We will print the results.

    print("The value of the constant beta is:")
    # We use sympy to format the output nicely
    beta_val = sympy.Rational(-1, 4)
    print(f"beta = {beta_val}")
    print("-" * 20)
    print("The expression for alpha^2 in terms of the constant curvature R is:")
    R = sympy.Symbol('R')
    kappa = sympy.Symbol('kappa')
    alpha_squared_expr = -R / (3 * kappa**2)

    # The problem asks to output each number in the final equation.
    # Let's format it to show the numerical coefficient clearly.
    coeff = alpha_squared_expr.coeff(R/kappa**2)

    print(f"alpha^2 = ({coeff}) * (R / kappa^2)")
    print("Which can be written as:")
    print(f"alpha^2 = -R / (3*kappa^2)")


if __name__ == "__main__":
    solve_supergravity_parameters()
    # Final answer submission format
    # Since the answer consists of two parts, a number for beta and an expression for alpha^2,
    # and the format requires a single value, we will present the numerical value for beta.
    final_answer_beta = -1/4
    # print(f'<<<{final_answer_beta}>>>')