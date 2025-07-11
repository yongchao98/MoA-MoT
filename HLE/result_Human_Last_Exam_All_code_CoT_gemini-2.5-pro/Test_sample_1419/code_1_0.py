import sympy

def find_wilson_fisher_fixed_point():
    """
    This script symbolically derives the leading order expression for the
    Wilson-Fisher fixed point coupling in phi^4 theory near d=4 dimensions.
    """

    # Define the symbols to be used in our equations.
    # u_star represents the fixed point coupling u*
    # u represents the generic coupling constant
    # epsilon is the small parameter in d = 4 - epsilon
    u_star, u, epsilon = sympy.symbols('u* u epsilon')
    pi = sympy.pi

    # --- Step 1 & 2: State the context and the beta function ---
    print("In ϕ⁴ theory near four dimensions (d = 4 - ε), the Renormalization Group (RG) flow of the coupling 'u' is described by the beta function.")
    print("To leading order (one-loop), the beta function is:")

    # Define the beta function expression
    beta_function_expr = -epsilon * u + (3 * u**2) / (16 * pi**2)
    beta_function_eq = sympy.Eq(sympy.Symbol('β(u)'), beta_function_expr)

    # Use sympy.pprint for a nicely formatted console output
    sympy.pprint(beta_function_eq, use_unicode=True)
    print("\n" + "="*50 + "\n")

    # --- Step 3: Explain the fixed point condition ---
    print(f"A fixed point, denoted as {u_star}, is a value of the coupling where the beta function is zero.")
    print(f"We must solve the equation β({u_star}) = 0:")

    # --- Step 4: Set up and solve the equation ---
    # Substitute u with u_star and set the expression to 0
    fixed_point_eq = sympy.Eq(beta_function_expr.subs(u, u_star), 0)
    sympy.pprint(fixed_point_eq, use_unicode=True)

    # Solve the equation for u_star
    solutions = sympy.solve(fixed_point_eq, u_star)

    print(f"\nThis equation has two solutions: {solutions}")

    # --- Step 5 & 6: Identify and present the non-trivial solution ---
    gaussian_fp = solutions[0]
    wilson_fisher_fp = solutions[1]

    print(f"The solution {gaussian_fp} corresponds to the trivial Gaussian fixed point (a free theory).")
    print("The non-trivial solution is the Wilson-Fisher fixed point, which describes critical phenomena.")

    # Extract the components of the final answer to print them explicitly
    coeff = wilson_fisher_fp / epsilon
    numerator = sympy.numer(coeff / pi**2)
    denominator = sympy.denom(coeff / pi**2)
    
    print("\nThe leading order expression for the Wilson-Fisher fixed point coupling is:")
    
    # Print the final equation with each number and symbol shown clearly
    final_eq_str = f"{u_star} = ({numerator} * π² / {denominator}) * ε"
    print("\n" + final_eq_str + "\n")


if __name__ == '__main__':
    find_wilson_fisher_fixed_point()