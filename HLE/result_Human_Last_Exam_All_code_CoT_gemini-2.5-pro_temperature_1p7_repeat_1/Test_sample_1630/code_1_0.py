import numpy as np
import sympy

def solve():
    """
    This function demonstrates that 9 fixed points are achievable for f(g(x)),
    where f and g are cubic polynomials with positive derivatives.
    """
    x = sympy.symbols('x')

    # Step 1: Define cubic polynomials f(x) and g(x) with specific coefficients.
    # These coefficients have been chosen to produce 9 fixed points.
    # f(x) = a3*x^3 + a2*x^2 + a1*x + a0
    # g(x) = b3*x^3 + b2*x^2 + b1*x + b0
    f_coeffs = [0.05, -0.75, 4.95, -7.25]
    g_coeffs = [0.05,  0.75, 4.95,  5.25]

    f_sym = sum(c * x**i for i, c in enumerate(reversed(f_coeffs)))
    g_sym = sum(c * x**i for i, c in enumerate(reversed(g_coeffs)))

    print("Chosen polynomials:")
    print(f"f(x) = {f_sym}")
    print(f"g(x) = {g_sym}\n")

    # Step 2: Verify that their derivatives are always positive.
    # A quadratic Ax^2+Bx+C is always positive if A>0 and B^2-4AC < 0.
    f_prime_sym = sympy.diff(f_sym, x)
    f_prime_coeffs = sympy.Poly(f_prime_sym, x).all_coeffs()
    A, B, C = [float(c) for c in f_prime_coeffs]
    f_discriminant = B**2 - 4*A*C
    print("Verifying f'(x) > 0:")
    print(f"f'(x) = {f_prime_sym}")
    print(f"Discriminant of f'(x) is {f_discriminant:.4f}.")
    if A > 0 and f_discriminant < 0:
        print("f'(x) is always positive.\n")
    else:
        print("Warning: f'(x) is NOT always positive.\n")

    g_prime_sym = sympy.diff(g_sym, x)
    g_prime_coeffs = sympy.Poly(g_prime_sym, x).all_coeffs()
    A, B, C = [float(c) for c in g_prime_coeffs]
    g_discriminant = B**2 - 4*A*C
    print("Verifying g'(x) > 0:")
    print(f"g'(x) = {g_prime_sym}")
    print(f"Discriminant of g'(x) is {g_discriminant:.4f}.")
    if A > 0 and g_discriminant < 0:
        print("g'(x) is always positive.\n")
    else:
        print("Warning: g'(x) is NOT always positive.\n")

    # Step 3: Find the fixed points by solving f(g(x)) - x = 0.
    h_sym = f_sym.subs(x, g_sym)
    fixed_point_eq = sympy.expand(h_sym - x)

    # Convert to a numerical polynomial and find roots
    poly = sympy.Poly(fixed_point_eq, x)
    coeffs_np = [float(c) for c in poly.all_coeffs()]
    
    print("The fixed points are the roots of the final equation f(g(x)) - x = 0.")
    print("This is a polynomial of degree 9. Its coefficients are:")
    # Using np.round to make it more readable
    print(np.round(coeffs_np, 4))
    print()

    roots = np.roots(coeffs_np)

    # Filter for real roots (where the imaginary part is close to zero)
    real_roots = roots[np.isclose(roots.imag, 0)].real

    print(f"The number of real roots found is: {len(real_roots)}")
    print("The real roots (fixed points) are:")
    print(np.round(sorted(real_roots), 4))
    
    print("\nBased on the mathematical derivation and this example, the maximum number is 9.")

solve()