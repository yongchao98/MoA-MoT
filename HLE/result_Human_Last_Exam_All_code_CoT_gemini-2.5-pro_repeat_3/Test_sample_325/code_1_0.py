import sympy

def solve_decoupling_exponent():
    """
    Calculates the sharp l^2 decoupling exponent for the curve
    gamma(t) = (cos(t), sin(t), t) in R^3.
    """

    # Step 1: Define the curve and the dimension of the space.
    n = 3
    t = sympy.Symbol('t')
    gamma = sympy.Matrix([sympy.cos(t), sympy.sin(t), t])

    # Step 2: Verify the non-degeneracy of the curve.
    # A curve in R^n is non-degenerate if its first n derivatives are linearly independent.
    # We check this by computing the determinant of the Wronskian matrix.
    gamma_p1 = sympy.diff(gamma, t)
    gamma_p2 = sympy.diff(gamma_p1, t)
    gamma_p3 = sympy.diff(gamma_p2, t)

    # Create the Wronskian matrix
    W = sympy.Matrix([gamma_p1.T, gamma_p2.T, gamma_p3.T]).T
    
    # Calculate the determinant
    det_W = W.det()
    # The determinant simplifies to cos(t)^2 + sin(t)^2 = 1.
    # Since the determinant is non-zero (it's 1), the curve is non-degenerate.
    
    # Step 3: State the theorem for the sharp exponent.
    # According to a theorem by Bourgain, Demeter, and Guth, the sharp l^2
    # decoupling exponent p_c for a non-degenerate curve in R^n (for n >= 2) is:
    # p_c = n(n+1) / (n-1)

    print("The sharp l^2 decoupling exponent p_c for a non-degenerate curve in R^n is given by the formula:")
    print("p_c = n * (n + 1) / (n - 1)\n")
    
    print(f"The given curve is gamma(t) = (cos(t), sin(t), t) in R^3, so the dimension n = {n}.")
    print("This curve is non-degenerate because the determinant of its first three derivatives is 1, which is non-zero.\n")

    # Step 4: Calculate the exponent for n=3.
    numerator_val = n * (n + 1)
    denominator_val = n - 1
    p_c = numerator_val / denominator_val
    
    print("Plugging n = 3 into the formula gives:")
    print(f"p_c = ({n} * ({n} + 1)) / ({n} - 1)")
    print(f"p_c = ({n} * {n + 1}) / {denominator_val}")
    print(f"p_c = {numerator_val} / {denominator_val}")
    print(f"p_c = {int(p_c)}\n")
    
    print(f"The sharp l^2 decoupling exponent is {int(p_c)}.")

solve_decoupling_exponent()