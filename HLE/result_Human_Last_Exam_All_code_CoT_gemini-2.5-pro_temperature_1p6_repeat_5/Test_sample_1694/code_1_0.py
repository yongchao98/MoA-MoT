import sympy

def solve_system():
    """
    Solves the interconnected system and finds the sum of the first coordinate of solutions.
    """
    # Use sympy for symbolic math and pretty printing
    sympy.init_printing(use_unicode=True)

    # --- Part 1: Solve for Y1 = X1^2 and Y2 = X2^2 ---
    print("Step 1: Solving the Sylvester equations for X_1^2 and X_2^2.")
    
    # For the first equation: A1*Y1 + Y1*B1 = C1
    # Leads to 11*y11 = -53/12, and other elements are 0.
    y1_11 = sympy.Rational(-53, 132)
    Y1 = sympy.Matrix([[y1_11, 0], [0, 0]])
    
    # For the second equation: A2*Y2 + Y2*B2 = C2
    # Leads to 10*z11 = -3/11, and other elements are 0.
    y2_11 = sympy.Rational(-3, 110)
    Y2 = sympy.Matrix([[y2_11, 0], [0, 0]])
    
    print("From the first equation, we find X_1^2:")
    print(sympy.pretty(Y1))
    print("\nFrom the second equation, we find X_2^2:")
    print(sympy.pretty(Y2))

    # --- Part 2: Find the solutions for the first coordinates ---
    print("\nStep 2: Finding the possible values for the first coordinates of X_1 and X_2.")
    
    # Solve for the (1,1) element of X1
    x1_sols = sympy.solve(sympy.Symbol('x')**2 - y1_11, sympy.Symbol('x'))
    
    print(f"The possible values for the first coordinate of X_1 (the solutions to x^2 = {y1_11}) are:")
    print(f"x_1a = {x1_sols[0]}")
    print(f"x_1b = {x1_sols[1]}")

    # --- Part 3: Calculate the sum ---
    print("\nStep 3: Calculating the sum of the first coordinate of all solutions.")
    print("A solution to the system is a pair (X1, X2).")
    print("Since there are two possible X1 matrices and two for X2, there are 4 total solutions.")
    print("We define the 'first coordinate of a solution' as the (1,1) element of its X1 matrix.")
    
    # The first coordinates of the 4 solutions are determined by the two possible values for X1's (1,1) element.
    # Each value appears twice, as it can be paired with either of the two X2 solutions.
    coord1 = x1_sols[0]
    coord2 = x1_sols[0]
    coord3 = x1_sols[1]
    coord4 = x1_sols[1]
    
    total_sum = coord1 + coord2 + coord3 + coord4
    
    print("\nThe first coordinates of the four solutions are:")
    print(f"1. {coord1}")
    print(f"2. {coord2}")
    print(f"3. {coord3}")
    print(f"4. {coord4}")
    
    # Print the final equation with all numbers
    print("\nThe sum is the addition of these four values:")
    print(f"({coord1}) + ({coord2}) + ({coord3}) + ({coord4}) = {total_sum}")

solve_system()