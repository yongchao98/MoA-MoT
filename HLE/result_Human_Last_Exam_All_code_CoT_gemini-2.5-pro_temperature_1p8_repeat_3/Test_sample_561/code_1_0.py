import math

def solve_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension of the surreal piano keys.
    """
    # 1. Define the parameters based on the problem description.
    N = 5  # Number of non-overlapping copies (5 black keys)

    # The transformation scales the 3x1 keyboard to fit into a black key.
    # Black key width = (1/2) * (3 units / 7 white keys) = 3/14
    # Black key height = 9/14 * (1 unit) = 9/14
    
    # The scaling factors (singular values) are the ratios of new dimensions to old.
    # We order them from largest to smallest for the formula.
    sigma_1 = (9/14) / 1  # Vertical scaling
    sigma_2 = (3/14) / 3  # Horizontal scaling

    # 2. Set up the dimension equation.
    # For a self-affine set with dimension 1 < D <= 2, the equation is:
    # N * sigma_1 * (sigma_2)^(D-1) = 1
    # This simplifies to N * sigma_1 = (1/sigma_2)^(D-1)
    # and further to N * sigma_1 / sigma_2 = (1/sigma_2)^D
    # which is equivalent to 45 = 14^D
    
    final_eq_lhs = 45
    final_eq_rhs_base = 14

    # 3. Solve the equation for D.
    # D = log_14(45) = log(45) / log(14)
    dimension = math.log(final_eq_lhs) / math.log(final_eq_rhs_base)

    print("The dimension D is found by solving the equation for self-affine fractals.")
    print("The general form of the equation is: N * sigma_1 * (sigma_2)^(D-1) = 1")
    print(f"With N={N}, sigma_1={sigma_1:.3f}, and sigma_2={sigma_2:.3f}, this simplifies to the final equation:")
    print(f"    {final_eq_lhs} = {final_eq_rhs_base}^D")
    print("\nSolving for D by taking the logarithm of both sides gives D = log(45) / log(14).")
    print("\nThe Minkowski–Bouligand dimension of the piano's black keys is:")
    print(dimension)

solve_fractal_dimension()