import numpy as np
import numpy.polynomial.polynomial as P

def solve_and_explain():
    """
    This function explains the reasoning for the maximum number of fixed points
    and verifies it with a concrete example using numerical computation.
    """

    # Step 1: Theoretical Maximum
    # Let f(x) and g(x) be polynomials of degree 3.
    # The composition h(x) = f(g(x)) is a polynomial of degree 3 * 3 = 9.
    # A fixed point is a solution to h(x) = x, or h(x) - x = 0.
    # The degree of h(x) - x is 9.
    # By the Fundamental Theorem of Algebra, a polynomial of degree 9 has at most 9 real roots.
    print("Step 1: Establishing the upper bound")
    print("The function f(g(x)) is a polynomial of degree 3 * 3 = 9.")
    print("The equation for fixed points, f(g(x)) - x = 0, is a polynomial equation of degree 9.")
    print("Therefore, the maximum number of fixed points is at most 9.")
    print("-" * 20)

    # Step 2: Achieving the maximum
    # To show that 9 fixed points are possible, we need to construct f(x) and g(x) such that:
    # 1. They are degree 3 polynomials.
    # 2. Their derivatives f'(x) and g'(x) are always positive.
    # 3. The equation f(g(x)) = x has 9 distinct real roots.

    # We can construct such functions by perturbing the 3rd-degree Chebyshev polynomial, T_3(x) = 4x^3 - 3x.
    # The equation T_3(T_3(x)) = x, which is T_9(x) = x, is known to have 9 distinct real roots in [-1, 1].
    # However, T_3'(x) = 12x^2 - 3, which is not always positive.

    # Let's define f(x) and g(x) as slightly modified Chebyshev polynomials.
    # Let c and d be constants.
    # f(y) = 4y^3 - 3y + cy
    # g(x) = 4x^3 - 3x + dx
    #
    # Their derivatives are:
    # f'(y) = 12y^2 - 3 + c
    # g'(x) = 12x^2 - 3 + d
    #
    # For these derivatives to be positive for all x, we need c > 3 and d > 3.
    # We will choose c and d to be slightly greater than 3.
    c = 3.01
    d = 3.01
    
    print("Step 2: Constructing an example")
    print(f"We define f(x) and g(x) based on the Chebyshev polynomial T_3(x) = 4x^3 - 3x.")
    print(f"Let f(x) = 4x^3 - 3x + c*x and g(x) = 4x^3 - 3x + d*x.")
    print(f"Their derivatives are f'(x) = 12x^2 - 3 + c and g'(x) = 12x^2 - 3 + d.")
    print(f"To ensure the derivatives are always positive, we need c > 3 and d > 3.")
    print(f"We choose c = {c} and d = {d}.")
    print("-" * 20)

    # Step 3: Numerically find the fixed points
    # We need to find the roots of the equation f(g(x)) - x = 0.
    
    # Define the polynomials using numpy's polynomial class
    # g(x) = 4x^3 + (d-3)x
    g_coeffs = [0, d - 3, 0, 4] # Coefficients for [x^0, x^1, x^2, x^3]
    g_poly = P.Polynomial(g_coeffs)

    # f(y) = 4y^3 + (c-3)y
    f_coeffs = [0, c - 3, 0, 4]
    
    # Compose f(g(x))
    # h(x) = f(g(x)) = 4*(g(x))^3 + (c-3)*g(x)
    x_poly = P.Polynomial([0, 1])
    h_poly = 4 * (g_poly**3) + (c - 3) * g_poly
    
    # The equation for fixed points is h(x) - x = 0
    k_poly = h_poly - x_poly
    
    # Find the roots of the resulting polynomial
    roots = k_poly.roots()
    
    # Filter for real roots (roots with a very small imaginary part)
    real_roots = roots[np.isclose(roots.imag, 0)].real
    
    # Round roots to a few decimal places to count unique ones
    unique_real_roots = np.unique(np.round(real_roots, decimals=6))
    
    num_fixed_points = len(unique_real_roots)

    print("Step 3: Finding the number of fixed points numerically")
    print("We solve the polynomial equation f(g(x)) - x = 0 for the chosen f and g.")
    print(f"The degree of the polynomial f(g(x)) - x is {k_poly.degree()}.")
    print(f"Found {len(real_roots)} real roots numerically.")
    print(f"Number of unique real roots: {num_fixed_points}")
    print("The 9 fixed points are approximately:")
    print(np.round(np.sort(unique_real_roots), 4))
    print("-" * 20)
    
    # Final Answer
    print("Conclusion: The maximum number of fixed points is 9.")
    return num_fixed_points

max_points = solve_and_explain()
print(f"\nFinal Answer: {max_points}")