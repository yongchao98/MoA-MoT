import sympy

def solve_and_sum():
    """
    Solves the two matrix equations and calculates the sum of the first coordinates of all solutions.
    """
    # Initialize symbols
    i = sympy.I

    # --- Equation 1 ---
    # A1*X1^2 + X1^2*B1 = C1
    # where A1 = diag(5, -5), B1 = diag(6, 6), C1 = diag(-53/12, 0)
    # Let Y1 = X1^2 = diag(y1, y2). The Sylvester equation gives:
    # (5+6)*y1 = -53/12  => 11*y1 = -53/12 => y1 = -53/132
    # (-5+6)*y2 = 0 => y2 = 0
    # So, X1^2 = diag(-53/132, 0)
    # Let X1 be a 2x2 matrix. The only solutions are diagonal with X1[1,1]=0 and X1[0,0]^2 = -53/132.
    y1_11 = sympy.Rational(-53, 132)
    x1_coords = sympy.sqrt(y1_11).as_real_imag()[1] * i, -sympy.sqrt(y1_11).as_real_imag()[1] * i

    # --- Equation 2 ---
    # A2*X2^2 + X2^2*B2 = C2
    # where A2 = diag(4, -5), B2 = diag(6, 6), C2 = diag(-3/11, 0)
    # Let Y2 = X2^2 = diag(y1, y2). The Sylvester equation gives:
    # (4+6)*y1 = -3/11 => 10*y1 = -3/11 => y1 = -3/110
    # (-5+6)*y2 = 0 => y2 = 0
    # So, X2^2 = diag(-3/110, 0)
    # Let X2 be a 2x2 matrix. The only solutions are diagonal with X2[1,1]=0 and X2[0,0]^2 = -3/110.
    y2_11 = sympy.Rational(-3, 110)
    x2_coords = sympy.sqrt(y2_11).as_real_imag()[1] * i, -sympy.sqrt(y2_11).as_real_imag()[1] * i
    
    # Combine all first coordinates
    all_coords = list(x1_coords) + list(x2_coords)
    
    # Calculate the sum
    total_sum = sum(all_coords)

    # Print the final equation as requested
    equation_str = " + ".join([f"({sympy.pretty(c, use_unicode=False)})" for c in all_coords])
    print("The solutions for the first coordinate of X_1 are:")
    print(f"x_1_a = {sympy.pretty(x1_coords[0], use_unicode=False)}")
    print(f"x_1_b = {sympy.pretty(x1_coords[1], use_unicode=False)}")
    print("\nThe solutions for the first coordinate of X_2 are:")
    print(f"x_2_a = {sympy.pretty(x2_coords[0], use_unicode=False)}")
    print(f"x_2_b = {sympy.pretty(x2_coords[1], use_unicode=False)}")
    print("\nThe sum of all these coordinates is:")
    print(f"{equation_str} = {sympy.pretty(total_sum, use_unicode=False)}")

solve_and_sum()