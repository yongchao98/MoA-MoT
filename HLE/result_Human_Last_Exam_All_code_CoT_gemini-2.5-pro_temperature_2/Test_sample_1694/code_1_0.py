import sympy

def solve_and_print():
    """
    This function solves the given system of matrix equations and
    finds the sum of the first coordinate of the solutions.
    """
    # Define the constant matrix B = 6I
    B = sympy.Matrix([[6, 0], [0, 6]])

    # --- Step 1 & 2: Solve the first equation for X1 ---
    print("--- Solving for X1 ---")
    
    # Define matrices from the first equation
    A1 = sympy.Matrix([[5, 0], [0, -5]])
    C1 = sympy.Matrix([[-sympy.Rational(53, 12), 0], [0, 0]])

    # The equation is (A1 + B) * X1^2 = C1
    # We solve for X1^2
    M1 = A1 + B
    X1_squared = M1.inv() * C1
    
    # The first element of X1^2 is x1_11^2
    y1_11 = X1_squared[0, 0]

    # As derived in the plan, the non-diagonal elements of X1 are zero.
    # The first coordinate of X1 is x1_11 = sqrt(y1_11).
    # This gives two solutions because of the square root.
    sol1_a = sympy.sqrt(y1_11)
    sol1_b = -sympy.sqrt(y1_11)

    print(f"The matrix X1^2 is calculated as:\n{sympy.pretty(X1_squared)}")
    print(f"The possible values for the first coordinate of X1 are:\n1: {sol1_a}\n2: {sol1_b}\n")

    # --- Step 3: Solve the second equation for X2 ---
    print("--- Solving for X2 ---")

    # Define matrices from the second equation
    A2 = sympy.Matrix([[4, 0], [0, -5]])
    C2 = sympy.Matrix([[-sympy.Rational(3, 11), 0], [0, 0]])
    
    # The equation is (A2 + B) * X2^2 = C2
    # We solve for X2^2
    M2 = A2 + B
    X2_squared = M2.inv() * C2
    
    # The first element of X2^2 is x2_11^2
    y2_11 = X2_squared[0, 0]

    # Similar to X1, the first coordinate of X2 has two solutions.
    sol2_a = sympy.sqrt(y2_11)
    sol2_b = -sympy.sqrt(y2_11)

    print(f"The matrix X2^2 is calculated as:\n{sympy.pretty(X2_squared)}")
    print(f"The possible values for the first coordinate of X2 are:\n1: {sol2_a}\n2: {sol2_b}\n")
    
    # --- Step 4: Calculate the total sum ---
    print("--- Calculating the Final Sum ---")
    
    # The total sum is the sum of all four first-coordinate solutions.
    total_sum = sol1_a + sol1_b + sol2_a + sol2_b

    # Print the sum as an explicit equation, showing each number.
    print("The final equation for the sum is:")
    print(f"({sol1_a}) + ({sol1_b}) + ({sol2_a}) + ({sol2_b}) = {total_sum}")

solve_and_print()