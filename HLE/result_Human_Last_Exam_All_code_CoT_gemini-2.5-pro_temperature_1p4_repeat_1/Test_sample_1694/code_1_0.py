import sympy

def solve_and_sum():
    """
    This function solves the given system of matrix equations and finds the sum
    of the first coordinate of all solution matrices.
    """
    # Initialize lists to store solution coordinates and the total sum
    all_first_coords = []
    total_sum = 0
    
    # --- Equation 1 ---
    print("--- Solving the first equation ---")
    
    # Define matrices for the first equation
    A1 = sympy.Matrix([[5, 0], [0, -5]])
    B1 = sympy.Matrix([[6, 0], [0, 6]])
    C1 = sympy.Matrix([[-sympy.Rational(53, 12), 0], [0, 0]])
    
    # Solve for Y1 = X1^2 in A1*Y1 + Y1*B1 = C1
    # Since all matrices are diagonal, we can solve for elements directly.
    # (A1_ii + B1_jj) * Y1_ij = C1_ij
    y1_11 = C1[0, 0] / (A1[0, 0] + B1[0, 0])
    y1_22 = C1[1, 1] / (A1[1, 1] + B1[1, 1])
    Y1 = sympy.Matrix([[y1_11, 0], [0, y1_22]])
    print(f"Calculated Y1 = X1^2:\n{sympy.pretty(Y1)}\n")

    # Solve for X1 by taking the square root of Y1.
    # The solutions for the first coordinate x1_11 are +/- sqrt(y1_11)
    x1_sol1 = sympy.sqrt(Y1[0, 0])
    x1_sol2 = -sympy.sqrt(Y1[0, 0])
    first_coords_x1 = [x1_sol1, x1_sol2]
    all_first_coords.extend(first_coords_x1)
    sum_x1 = sum(first_coords_x1)

    print("The two solutions for the first coordinate of X1 are:")
    print(f"1. {sympy.pretty(x1_sol1)}")
    print(f"2. {sympy.pretty(x1_sol2)}\n")

    # --- Equation 2 ---
    print("--- Solving the second equation ---")
    
    # Define matrices for the second equation
    A2 = sympy.Matrix([[4, 0], [0, -5]])
    B2 = sympy.Matrix([[6, 0], [0, 6]])
    C2 = sympy.Matrix([[-sympy.Rational(3, 11), 0], [0, 0]])
    
    # Solve for Y2 = X2^2 in A2*Y2 + Y2*B2 = C2
    y2_11 = C2[0, 0] / (A2[0, 0] + B2[0, 0])
    y2_22 = C2[1, 1] / (A2[1, 1] + B2[1, 1])
    Y2 = sympy.Matrix([[y2_11, 0], [0, y2_22]])
    print(f"Calculated Y2 = X2^2:\n{sympy.pretty(Y2)}\n")

    # Solve for X2 by taking the square root of Y2.
    # The solutions for the first coordinate x2_11 are +/- sqrt(y2_11)
    x2_sol1 = sympy.sqrt(Y2[0, 0])
    x2_sol2 = -sympy.sqrt(Y2[0, 0])
    first_coords_x2 = [x2_sol1, x2_sol2]
    all_first_coords.extend(first_coords_x2)
    sum_x2 = sum(first_coords_x2)
    
    print("The two solutions for the first coordinate of X2 are:")
    print(f"1. {sympy.pretty(x2_sol1)}")
    print(f"2. {sympy.pretty(x2_sol2)}\n")

    # --- Final Summation ---
    print("--- Calculating the total sum ---")
    
    total_sum = sum(all_first_coords)
    
    # Format the final equation string
    equation_parts = []
    for coord in all_first_coords:
        equation_parts.append(f"({sympy.pretty(coord).replace('/n', '')})")
    equation_str = " + ".join(equation_parts)
    
    print("The sum of the first coordinates of all solutions is given by the equation:")
    print(f"{equation_str} = {total_sum}")

if __name__ == "__main__":
    solve_and_sum()
<<<0>>>