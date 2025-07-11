import numpy as np

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # Step 1: Solve the homogeneous equation to find B and D
    # Characteristic equation: 8r^2 - 6r + 1 = 0
    char_eq_coeffs = [8, -6, 1]
    roots = np.roots(char_eq_coeffs)

    # Assign B to the larger root and D to the smaller root, as is conventional.
    B = max(roots)
    D = min(roots)

    # Step 2: Find the particular solution E
    # 8E - 6E + E = 1 => 3E = 1
    E = 1/3

    # Step 3: Use initial conditions to find A and C
    y0 = 1.0
    ym1 = 2.0

    # Calculate y[1] from the difference equation: 8y[1] - 6y[0] + y[-1] = 1
    # 8*y[1] - 6*(1) + 2 = 1  => 8*y[1] - 4 = 1 => 8*y[1] = 5
    y1 = 5.0 / 8.0

    # Set up the system of linear equations for A and C
    # Equation 1 (from y[0]): A + C = y[0] - E
    # Equation 2 (from y[1]): A*B + C*D = y[1] - E
    matrix = np.array([[1, 1], [B, D]])
    vector = np.array([y0 - E, y1 - E])

    # Solve for A and C
    try:
        A, C = np.linalg.solve(matrix, vector)
    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")
        return

    # Step 4: Print the parameters and calculate the final expression
    print("The closed form solution is of the form: y[n] = A * (B)^n + C * (D)^n + E")
    print("The determined parameters are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print("\n")
    print("The final equation is:")
    print(f"y[n] = {A} * ({B})^n + {C} * ({D})^n + {E}")
    print("\n")

    result = E / A + (D * C) / B
    print(f"The value of E/A + (D*C)/B is: {result}")


if __name__ == "__main__":
    solve_difference_equation()
<<<0.75>>>