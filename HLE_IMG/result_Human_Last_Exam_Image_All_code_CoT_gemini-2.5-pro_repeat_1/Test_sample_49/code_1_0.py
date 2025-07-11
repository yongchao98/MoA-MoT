import sympy

def solve_cutoff_frequency():
    """
    This function symbolically calculates the cutoff frequency of the given ladder network.
    """
    # Define the symbols for resistance, capacitance, and the ladder's characteristic resistance
    r, C, R_L = sympy.symbols('r C R_L', positive=True)

    # Step 1: Set up and solve the equation for the characteristic resistance R_L
    # The equation is R_L = r || (2r + R_L)
    # R_L = (r * (2*r + R_L)) / (r + (2*r + R_L))
    # R_L * (3*r + R_L) = 2*r**2 + r*R_L
    # R_L**2 + 2*r*R_L - 2*r**2 = 0
    equation = sympy.Eq(R_L**2 + 2*r*R_L - 2*r**2, 0)

    # Solve the quadratic equation for R_L
    R_L_solutions = sympy.solve(equation, R_L)

    # The resistance must be positive, so we select the positive solution
    R_L_expr = [sol for sol in R_L_solutions if sol.is_positive][0]
    
    print("Step 1: Calculating the characteristic resistance of the infinite ladder (R_L)")
    print(f"The self-consistency equation is: R_L^2 + 2*r*R_L - 2*r^2 = 0")
    print(f"The positive solution for R_L is: {R_L_expr}\n")

    # Step 2: Calculate the Thevenin equivalent resistance R_th
    # R_th = r + R_L + r = 2r + R_L
    R_th_expr = 2*r + R_L_expr
    R_th_simplified = sympy.simplify(R_th_expr)

    print("Step 2: Calculating the Thevenin resistance (R_th) seen by the capacitor")
    print(f"R_th = 2*r + R_L")
    print(f"Substituting R_L, we get R_th = {R_th_simplified}\n")

    # Step 3: Determine the cutoff frequency f_c
    # f_c = 1 / (2 * pi * R_th * C)
    pi = sympy.pi
    f_c_expr = 1 / (2 * pi * R_th_simplified * C)
    
    print("Step 3: Determining the cutoff frequency (f_c)")
    print("f_c = 1 / (2 * pi * R_th * C)")
    print(f"The final expression for the cutoff frequency is:")
    # The default sympy printing can be a bit messy, let's format it.
    # R_th = r*(1 + sqrt(3))
    # f_c = 1 / (2 * pi * C * r * (1 + sqrt(3)))
    print(f"f_c = 1 / (2 * \u03C0 * C * r * (1 + \u221A3))")
    
    print("\nThe numbers in the final equation f_c = 1 / (2 * \u03C0 * C * r * (A + \u221AB)) are:")
    # Extracting coefficients from sympy expression can be complex. We know them from derivation.
    # R_th_simplified = r * (sqrt(3) + 1)
    # So A=1, B=3
    print("Numerator: 1")
    print("Denominator term (2 * pi * ...): 2")
    print("A = 1")
    print("B = 3")

solve_cutoff_frequency()