import numpy as np

def solve_and_sum_coordinates():
    """
    Solves the given system of matrix equations and calculates the sum of the
    first coordinates of the solutions.
    """
    # === Equation 1 ===
    # M1 * X1^2 + X1^2 * N1 = R1
    M1 = np.array([[5, 0], [0, -5]])
    N1 = np.array([[6, 0], [0, 6]])
    R1 = np.array([[-53/12, 0], [0, 0]])

    # The equation simplifies to (M1 + N1) * X1^2 = R1
    # Let Y1 = X1^2. Then Y1 = inv(M1 + N1) * R1
    M1_plus_N1 = M1 + N1
    inv_M1_plus_N1 = np.linalg.inv(M1_plus_N1)
    Y1 = inv_M1_plus_N1 @ R1
    
    # === Equation 2 ===
    # M2 * X2^2 + X2^2 * N2 = R2
    M2 = np.array([[4, 0], [0, -5]])
    N2 = np.array([[6, 0], [0, 6]]) # N2 is the same as N1
    R2 = np.array([[-3/11, 0], [0, 0]])
    
    # The equation simplifies to (M2 + N2) * X2^2 = R2
    # Let Y2 = X2^2. Then Y2 = inv(M2 + N2) * R2
    M2_plus_N2 = M2 + N2
    inv_M2_plus_N2 = np.linalg.inv(M2_plus_N2)
    Y2 = inv_M2_plus_N2 @ R2

    # Y1 and Y2 are diagonal. We need to find their matrix square roots.
    # For a diagonal matrix D = diag(d1, d2), a square root is S = diag(s1, s2)
    # where s1^2 = d1 and s2^2 = d2.
    
    # The (1,1) element of Y1 is Y1[0, 0]. The square roots of this element will be the
    # (1,1) elements of the solution matrices X1.
    y1_11 = Y1[0, 0]
    
    # np.sqrt handles complex numbers correctly.
    x1_11_sol1 = np.sqrt(y1_11)
    x1_11_sol2 = -x1_11_sol1

    # For X2, the solutions are independent.
    # The problem has 2 solutions for X1 and 2 for X2, leading to 4 total solution pairs (X1, X2).
    # Solution 1: (X1_sol1, X2_sol1)
    # Solution 2: (X1_sol1, X2_sol2)
    # Solution 3: (X1_sol2, X2_sol1)
    # Solution 4: (X1_sol2, X2_sol2)
    
    # The "first coordinate of a solution" is interpreted as the (1,1) element of the X1 matrix.
    # For the 4 solutions, these coordinates are:
    coords = [x1_11_sol1, x1_11_sol1, x1_11_sol2, x1_11_sol2]
    
    # Summing these coordinates
    total_sum = sum(coords)

    print("The system has 4 solutions. The first coordinate of each solution (the (1,1) entry of X1) is determined by the two possible square roots of the (1,1) entry of X1^2.")
    print(f"The (1,1) entry of X1^2 is {y1_11:.4f}.")
    print(f"The two possible values for the first coordinate of X1 are {x1_11_sol1} and {x1_11_sol2}.")
    print("\nThe four solution pairs (X1, X2) result in the following four values for the first coordinate:")
    
    coord_strings = [f"({c.real:.4f}{c.imag:+.4f}j)" for c in coords]
    print(f"{coord_strings[0]}, {coord_strings[1]}, {coord_strings[2]}, {coord_strings[3]}")
    
    print("\nThe sum of the first coordinate of solutions is:")
    sum_string = f" {coord_strings[0]} + {coord_strings[1]} + {coord_strings[2]} + {coord_strings[3]}"
    print(f"{sum_string} = {total_sum}")

    # The result is 0, so we print its real part.
    print(f"\nFinal Sum: {total_sum.real}")

solve_and_sum_coordinates()