import sympy

def find_wilson_fisher_fixed_point():
    """
    This function symbolically derives the leading order expression for the
    Wilson-Fisher fixed point in phi^4 theory near d=4 dimensions.
    """
    # 1. Define the necessary mathematical symbols
    # u: the coupling constant
    # epsilon: the small parameter in d = 4 - epsilon
    # u_star: the symbol for the fixed point coupling for printing
    u = sympy.Symbol('u')
    epsilon = sympy.Symbol('ϵ')
    u_star = sympy.Symbol('u^*')
    pi = sympy.pi

    # 2. Define the one-loop beta function for the phi^4 theory
    # beta(u) = -epsilon*u + 3*u^2 / (8*pi^2)
    beta_function = -epsilon * u + (3 * u**2) / (8 * pi**2)

    # 3. Solve for the fixed points where beta(u) = 0
    fixed_points = sympy.solve(beta_function, u)

    # 4. Identify the non-trivial Wilson-Fisher fixed point (the one not equal to zero)
    wilson_fisher_fp = None
    for fp in fixed_points:
        if fp != 0:
            wilson_fisher_fp = fp
            break

    # 5. Print the final result as a formatted equation,
    # highlighting the numbers involved as requested.
    if wilson_fisher_fp is not None:
        # The expression is (8/3)*pi**2*epsilon
        eight = 8
        three = 3
        two = 2
        print("The leading order expression for the fixed point coupling is:")
        print(f"{u_star} = ({eight} * {pi}**{two} / {three}) * {epsilon}")
    else:
        print("Could not find the non-trivial fixed point.")

if __name__ == "__main__":
    find_wilson_fisher_fixed_point()
