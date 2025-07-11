import sympy

def solve_matrix_system():
    """
    Solves the given system of matrix equations to find the sum of the
    first coordinates of all solution matrices.
    """
    # --- Step 1: Solve the first equation for X1 ---

    # The first equation is:
    # A1 * X1^2 + X1^2 * B = C1
    # Let Y1 = X1^2. The equation is A1*Y1 + Y1*B = C1.
    # Since B = 6*I, this simplifies to (A1 + 6*I) * Y1 = C1.
    # The matrices are diagonal, so we can solve for the (1,1) element of Y1.
    
    # From the problem statement for the first equation:
    a1_11 = sympy.Integer(5)
    c1_11 = sympy.Rational(-53, 12)
    b_val = sympy.Integer(6)
    
    # Calculate the (1,1) element of Y1 = X1^2
    y1_11 = c1_11 / (a1_11 + b_val)
    
    # We found X1^2 must be of the form [[y1_11, 0], [0, 0]].
    # Taking the square root implies that for X1 = [[a, b], [c, d]],
    # we must have b=c=d=0 and a^2 = y1_11.
    
    # Find the possible values for 'a', the first coordinate of X1
    a = sympy.Symbol('a')
    solutions_x1 = sympy.solve(a**2 - y1_11, a)
    sol1_a1 = solutions_x1[0]
    sol1_a2 = solutions_x1[1]

    # --- Step 2: Solve the second equation for X2 ---
    
    # The second equation has the same structure:
    # (A2 + 6*I) * Y2 = C2
    
    # From the problem statement for the second equation:
    a2_11 = sympy.Integer(4)
    c2_11 = sympy.Rational(-3, 11)

    # Calculate the (1,1) element of Y2 = X2^2
    y2_11 = c2_11 / (a2_11 + b_val)
    
    # Similarly, for X2 = [[a', b'], [c', d']], we find
    # b'=c'=d'=0 and a'^2 = y2_11.
    
    # Find the possible values for 'a_prime', the first coordinate of X2
    a_prime = sympy.Symbol('a_prime')
    solutions_x2 = sympy.solve(a_prime**2 - y2_11, a_prime)
    sol2_a1 = solutions_x2[0]
    sol2_a2 = solutions_x2[1]

    # --- Step 3: Sum the coordinates and print the final equation ---
    
    # The total sum is the sum of all found first coordinates
    total_sum = sol1_a1 + sol1_a2 + sol2_a1 + sol2_a2
    
    print("The solutions for the first coordinate of X1 are:")
    print(f"  {sol1_a1} and {sol1_a2}")
    
    print("\nThe solutions for the first coordinate of X2 are:")
    print(f"  {sol2_a1} and {sol2_a2}")

    print("\nThe final equation for the sum of all first coordinates is:")
    # The f-string formatting will represent the complex numbers appropriately.
    print(f"({sol1_a1}) + ({sol1_a2}) + ({sol2_a1}) + ({sol2_a2}) = {total_sum}")

# Execute the solver function
solve_matrix_system()