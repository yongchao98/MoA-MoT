import sympy as sp
import numpy as np

def solve_quadric_intersections():
    """
    Finds the number of real intersection points of two quadrics by using a
    coordinate transformation to simplify the system and then solving it algebraically.
    """
    # Define symbolic variables
    x, y, u, v = sp.symbols('x y u v')

    # Define the two quadric equations from the user's prompt
    Q1_expr = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    Q2_expr = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    # Define the inverse coordinate transformation: x = u, y = (v + 3u) / 2
    transform = {x: u, y: (v + 3*u) / 2}

    # Apply the transformation to the first equation and simplify
    # This simplifies to 2*u**2 - 16*u + 18*v**2 + 31, which can be written as:
    Q1_uv_std = 2*(u - 4)**2 + 18*v**2 - 1

    # Apply the transformation to the second equation and simplify
    Q2_uv = sp.simplify(Q2_expr.subs(transform))

    # From the simplified Q1, we can express 18*v**2
    # 18*v**2 = 1 - 2*(u - 4)**2
    # We can use this to establish a simpler relationship between u and v from Q2
    # The transformed Q2 is: 9*u**2 - 42*u*v - 74*u + 81*v**2 + 162*v + 149 = 0
    # Substitute 81*v**2 = (9/2) * (1 - 2*(u-4)**2) into it.
    # This leads to: -42*u*v - 2*u + 162*v + 9.5 = 0
    # Solving for v gives:
    v_in_terms_of_u = (2*u - 9.5) / (162 - 42*u)

    # Substitute this expression for v back into the simplified ellipse equation
    final_equation_in_u = Q1_uv_std.subs(v, v_in_terms_of_u)

    # To solve final_equation_in_u = 0, we find the roots of its numerator
    numerator_poly, denominator_poly = sp.fraction(sp.simplify(final_equation_in_u))

    # To ensure integer coefficients for printing, we can clear the denominators (from 9.5)
    # The expression for the numerator is:
    # 72*(u-4)**2*(27-7*u)**2 + 18*(2*u-9.5)**2 - 36*(27-7*u)**2
    # Multiplying by 4 clears the decimal from 9.5
    final_poly_to_solve = sp.expand(numerator_poly * 4)

    # Find the real roots of this polynomial
    poly_coeffs = [float(c) for c in sp.Poly(final_poly_to_solve, u).all_coeffs()]
    all_roots = np.roots(poly_coeffs)
    real_roots_u = [r.real for r in all_roots if abs(r.imag) < 1e-9]

    # For intersection points to be real, u must be in a specific range
    # derived from 2*(u - 4)**2 <= 1
    u_min = 4 - 1/np.sqrt(2)
    u_max = 4 + 1/np.sqrt(2)

    # Filter the roots to find those that correspond to real intersection points
    valid_roots = sorted([r for r in real_roots_u if u_min <= r <= u_max])
    
    # The number of distinct valid roots gives the number of intersection points
    # We use a tolerance to find distinct roots
    distinct_valid_roots = []
    if valid_roots:
        distinct_valid_roots.append(valid_roots[0])
        for i in range(1, len(valid_roots)):
            if abs(valid_roots[i] - valid_roots[i-1]) > 1e-9:
                distinct_valid_roots.append(valid_roots[i])

    num_intersections = len(distinct_valid_roots)

    # Print the explanation and the final polynomial equation
    print("The system of two quadric equations is solved by a change of variables.")
    print("This simplifies the problem to finding the roots of a single polynomial.")
    print("\nThe final polynomial equation in the variable u is:")
    
    # Use sympy's pretty print to display the equation
    final_poly_str = sp.pretty(sp.Poly(final_poly_to_solve, u), use_unicode=True)
    print(f"{final_poly_str} = 0")
    
    print(f"\nThis is a 4th-degree polynomial. Its real roots are approximately:")
    print([round(r, 4) for r in real_roots_u])
    
    print(f"\nFor these roots to correspond to real intersection points, they must be in the interval [{u_min:.4f}, {u_max:.4f}].")
    print(f"The distinct valid roots in this interval are approximately:")
    print([round(r, 4) for r in distinct_valid_roots])
    
    print(f"\nThere are {num_intersections} distinct valid roots, which means there are {num_intersections} real intersection points.")

solve_quadric_intersections()