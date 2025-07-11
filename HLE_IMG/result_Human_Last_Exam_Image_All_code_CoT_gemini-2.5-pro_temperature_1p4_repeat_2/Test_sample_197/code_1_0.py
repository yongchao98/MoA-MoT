import sympy as sp

def solve_and_explain():
    """
    This function performs the symbolic math to identify the correct curve.
    """
    # Define the variable x
    x = sp.Symbol('x')

    # Step 1: Define the function f(x) based on graphical analysis.
    # From the graph, we deduce:
    # - Vertical asymptote at x=2
    # - Local max at (1, 0) and local min at (3, 6)
    # This leads to the function f(x) = 1.5*x + 1.5/(x-2)
    a = 1.5
    b = 0
    c = 1.5
    d = 2
    f = a*x + b + c/(x - d)
    
    # Step 2: Calculate the second derivative of f(x).
    f_prime = sp.diff(f, x)
    f_double_prime = sp.diff(f_prime, x)

    # Step 3: Construct the target function g(x) = -0.5f''(3x-2)+1.
    # The equation to find is y = -0.5*f''(3x-2)+1
    # First, substitute (3x-2) into f''(x).
    f_double_prime_transformed = f_double_prime.subs(x, 3*x - 2)
    
    # Then, construct the final function g(x).
    scale = -0.5
    shift = 1
    g = scale * f_double_prime_transformed + shift

    # Step 4: Analyze the properties of g(x).
    # Find the horizontal asymptote by taking the limit as x -> oo.
    horizontal_asymptote = sp.limit(g, x, sp.oo)

    # Find the vertical asymptote by finding the root of the denominator.
    # The denominator is part of f_double_prime_transformed.
    # We find where the function is undefined.
    denominator = sp.fraction(f_double_prime_transformed)[1]
    vertical_asymptote_eq = sp.solve(denominator, x)

    # Step 5: Print the results and the conclusion.
    print("The function to identify is: y = -0.5 * f''(3*x - 2) + 1\n")
    print(f"Step 1: From the graph, we deduce the blue function is f(x) = {sp.pretty(f)}.\n")
    print(f"Step 2: The second derivative is f''(x) = {sp.pretty(f_double_prime)}.\n")
    print("Step 3: We construct the target function.")
    print(f"Substituting (3*x - 2) into f''(x) gives: f''(3x - 2) = {sp.pretty(f_double_prime_transformed)}")
    print(f"So, y = {scale} * ({sp.pretty(f_double_prime_transformed)}) + {shift}")
    final_eq = sp.simplify(g)
    print(f"The final simplified equation is: y = {sp.pretty(final_eq)}\n")

    print("Step 4: Analyze the properties of the resulting function.")
    print(f"Horizontal Asymptote: y = {horizontal_asymptote}")
    print(f"Vertical Asymptote: x = {vertical_asymptote_eq[0]} ≈ {float(vertical_asymptote_eq[0]):.2f}\n")

    print("Step 5: Match properties to the curves.")
    print("The calculated function has a horizontal asymptote at y=1 and a vertical asymptote at x=4/3 (~1.33).")
    print("Observing the graph:")
    print("- The red curve has a horizontal asymptote at y = -1.")
    print("- The green curve has a horizontal asymptote at y = 1, but its vertical asymptote is at x < 1.")
    print("- The purple curve has a horizontal asymptote at y = 2.")
    print("- The black curve has a horizontal asymptote at y = 1 and a vertical asymptote between x=1 and x=2.")
    print("\nConclusion: The black curve matches these properties.")

solve_and_explain()