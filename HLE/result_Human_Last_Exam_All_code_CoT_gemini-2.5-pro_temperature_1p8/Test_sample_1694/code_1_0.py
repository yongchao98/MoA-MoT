import sympy

def solve_and_sum():
    """
    Solves the given system of matrix equations and finds the sum of the first coordinates of the solutions.
    """
    # Step 1: Define matrices for the first equation
    A1 = sympy.Matrix([[5, 0], [0, -5]])
    B = sympy.Matrix([[6, 0], [0, 6]])
    C1 = sympy.Matrix([[-sympy.Rational(53, 12), 0], [0, 0]])
    I = sympy.eye(2)

    # Solve for Y1 = X1**2 from (A1 + 6*I)*Y1 = C1
    LHS1_matrix = A1 + 6 * I
    Y1 = LHS1_matrix.inv() * C1
    y1_11 = Y1[0, 0]

    # Find the possible values for the first coordinate of X1
    x1_sol_1 = sympy.sqrt(y1_11)
    x1_sol_2 = -sympy.sqrt(y1_11)

    print("Step 1: Solve for the first coordinate of X1.")
    print(f"The first matrix equation leads to X1^2 = {Y1.tolist()}.")
    print(f"The possible values for the first coordinate of X1 are:\nx1_a = {x1_sol_1}\nx1_b = {x1_sol_2}\n")

    # Step 2: Define matrices for the second equation
    A2 = sympy.Matrix([[4, 0], [0, -5]])
    C2 = sympy.Matrix([[-sympy.Rational(3, 11), 0], [0, 0]])

    # Solve for Y2 = X2**2 from (A2 + 6*I)*Y2 = C2
    LHS2_matrix = A2 + 6 * I
    Y2 = LHS2_matrix.inv() * C2
    y2_11 = Y2[0, 0]

    # Find the possible values for the first coordinate of X2
    x2_sol_1 = sympy.sqrt(y2_11)
    x2_sol_2 = -sympy.sqrt(y2_11)

    print("Step 2: Solve for the first coordinate of X2.")
    print(f"The second matrix equation leads to X2^2 = {Y2.tolist()}.")
    print(f"The possible values for the first coordinate of X2 are:\nx2_a = {x2_sol_1}\nx2_b = {x2_sol_2}\n")

    # Step 3: Sum the first coordinates over all four solutions
    # The four solution pairs for (x1_11, x2_11) give four terms to sum
    term1 = x1_sol_1 + x2_sol_1
    term2 = x1_sol_1 + x2_sol_2
    term3 = x1_sol_2 + x2_sol_1
    term4 = x1_sol_2 + x2_sol_2

    total_sum = term1 + term2 + term3 + term4

    print("Step 3: Calculate the final sum.")
    print("We sum the value (x1_11 + x2_11) for all four solution pairs.")
    print("The terms of the final sum are:")
    print(f"Term 1: ({x1_sol_1}) + ({x2_sol_1}) = {term1}")
    print(f"Term 2: ({x1_sol_1}) + ({x2_sol_2}) = {term2}")
    print(f"Term 3: ({x1_sol_2}) + ({x2_sol_1}) = {term3}")
    print(f"Term 4: ({x1_sol_2}) + ({x2_sol_2}) = {term4}")
    print("\nThe final equation for the sum is:")
    print(f"Sum = ({term1}) + ({term2}) + ({term3}) + ({term4})")
    print(f"Sum = {total_sum}")


if __name__ == "__main__":
    solve_and_sum()