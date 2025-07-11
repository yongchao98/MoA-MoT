import cmath

def solve_and_sum():
    """
    Solves the two matrix equations for the first coordinate of the solutions
    and calculates their sum.
    """
    # Equation 1: (A1 + 6I) * X1^2 = C1
    # We only need the (1,1) element of each matrix.
    a1_11 = 5
    b_11 = 6
    c1_11 = -53/12

    # (a1_11 + b_11) * x1_11^2 = c1_11
    # 11 * x1_11^2 = -53/12
    x1_11_squared = c1_11 / (a1_11 + b_11)

    # Find the two solutions for x1_11
    sol1_x1 = cmath.sqrt(x1_11_squared)
    sol2_x1 = -sol1_x1

    # Equation 2: (A2 + 6I) * X2^2 = C2
    # We only need the (1,1) element of each matrix.
    a2_11 = 4
    # b_11 is the same (6)
    c2_11 = -3/11

    # (a2_11 + b_11) * x2_11^2 = c2_11
    # 10 * x2_11^2 = -3/11
    x2_11_squared = c2_11 / (a2_11 + b_11)

    # Find the two solutions for x2_11
    sol1_x2 = cmath.sqrt(x2_11_squared)
    sol2_x2 = -sol1_x2

    # The four solutions for the first coordinate are sol1_x1, sol2_x1, sol1_x2, sol2_x2.
    total_sum = sol1_x1 + sol2_x1 + sol1_x2 + sol2_x2

    # Print the results as requested
    print("The four solutions for the first coordinates are:")
    print(f"s1 = {sol1_x1}")
    print(f"s2 = {sol2_x1}")
    print(f"s3 = {sol1_x2}")
    print(f"s4 = {sol2_x2}")
    print("\nThe final equation for the sum is:")
    # Using {:.4f} for cleaner display of complex numbers
    print(f"({sol1_x1:.4f}) + ({sol2_x1:.4f}) + ({sol1_x2:.4f}) + ({sol2_x2:.4f}) = {total_sum}")

solve_and_sum()