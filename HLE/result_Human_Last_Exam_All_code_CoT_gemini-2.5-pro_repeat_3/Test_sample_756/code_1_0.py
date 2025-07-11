import sympy

def solve():
    """
    This function finds the maximum value of |b| + |c| based on the derivation.
    """
    # Step 1: Define the parameter xv (the vertex position) and the expression for |b|+|c|.
    # From the derivation, for a family of valid polynomials, |b|+|c| can be expressed
    # as a function of xv.
    xv = sympy.Symbol('xv')
    f = (-xv**2 + 6*xv + 1) / (1 + xv)**2

    # Step 2: Find the maximum of this function using calculus.
    # Take the derivative with respect to xv.
    f_prime = sympy.diff(f, xv)

    # Solve for the critical points by setting the derivative to zero.
    critical_points = sympy.solve(f_prime, xv)

    # We select the critical point that is in our valid interval [0, 1).
    optimal_xv = None
    for p in critical_points:
        # Check if the point is a number and is in the interval [0, 1)
        if p.is_real and 0 <= p < 1:
            optimal_xv = p
            break
    
    if optimal_xv is None:
        print("Could not find an optimal point in the specified interval.")
        return

    # Step 3: For the optimal xv, find the specific coefficients b and c.
    # b = 4*xv / (1 + xv)**2
    # c = (1 + 2*xv - xv**2) / (1 + xv)**2
    b = (4*optimal_xv) / (1 + optimal_xv)**2
    c = (1 + 2*optimal_xv - optimal_xv**2) / (1 + optimal_xv)**2

    # Step 4: Calculate the final maximum value and print the details of the calculation.
    abs_b = abs(b)
    abs_c = abs(c)
    result = abs_b + abs_c
    
    a = -2 / (1 + optimal_xv)**2
    
    print(f"The maximum value is achieved for the polynomial with vertex at xv = {optimal_xv}.")
    print(f"The coefficients of this polynomial P(x) = ax^2 + bx + c are:")
    print(f"a = {a} (which is {float(a):.4f})")
    print(f"b = {b} (which is {float(b):.4f})")
    print(f"c = {c} (which is {float(c):.4f})")
    print("\nThe calculation for the maximum value of |b| + |c| is as follows:")
    print(f"|b| + |c| = |{b}| + |{c}|")
    print(f"          = {abs_b} + {abs_c}")
    print(f"          = {result}")

solve()