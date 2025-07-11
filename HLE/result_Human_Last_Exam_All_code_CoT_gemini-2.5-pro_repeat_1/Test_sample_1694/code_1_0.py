import numpy as np

def solve_and_sum_first_coordinates():
    """
    Solves the matrix system and calculates the sum of the first coordinates of the solutions.
    """
    # The system of equations is:
    # 1) (5 0; 0 -5) X1^2 + X1^2 (6 0; 0 6) = (-53/12 0; 0 0)
    # 2) (4 0; 0 -5) X2^2 + X2^2 (6 0; 0 6) = (-3/11 0; 0 0)
    # The two equations are independent.

    # Let's solve the first equation for X1.
    # A * X1^2 + X1^2 * B = C1
    # Since B = 6I, this simplifies to (A + 6I) * X1^2 = C1
    A = np.array([[5, 0], [0, -5]])
    C1 = np.array([[-53/12, 0], [0, 0]], dtype=float)
    I = np.identity(2)

    # Calculate K1 = A + 6I
    K1 = A + 6 * I

    # Solve for X1^2 = inv(K1) * C1
    try:
        inv_K1 = np.linalg.inv(K1)
        X1_squared = inv_K1 @ C1
    except np.linalg.LinAlgError:
        print("Matrix K1 is singular, cannot solve for X1^2.")
        return

    # From our analysis, X1 must be a matrix whose square is X1_squared.
    # X1^2 = [[x^2, 0], [0, 0]]
    # The top-left element of X1_squared gives us the square of the first coordinate of X1.
    lambda1 = X1_squared[0, 0]

    # The first coordinate of the X1 solutions are the square roots of lambda1.
    # Since lambda1 is negative, the coordinates will be imaginary.
    x1_coord1 = np.sqrt(np.complex(lambda1))
    x1_coord2 = -x1_coord1

    # The second equation for X2 can be solved independently and will also yield two solutions,
    # X2 and -X2, as it has the same structure.
    # This means for each of the 2 solutions for X1, there are 2 solutions for X2.
    # Total solutions for the system (X1, X2) is 2 * 2 = 4.
    # The solutions are (X1_sol1, X2_sol1), (X1_sol1, X2_sol2), (X1_sol2, X2_sol1), (X1_sol2, X2_sol2).

    # We are asked to find the sum of the first coordinate of solutions.
    # We interpret this as the sum of the first coordinate of X1 for each of the 4 solutions.
    # The first coordinates of X1 for these four solutions are:
    coords = [x1_coord1, x1_coord1, x1_coord2, x1_coord2]

    # Calculate the sum
    total_sum = sum(coords)

    # Print the final equation with each number
    print("The first coordinates of X1 for the four solutions of the system are:")
    print(f"Solution 1: {coords[0]}")
    print(f"Solution 2: {coords[1]}")
    print(f"Solution 3: {coords[2]}")
    print(f"Solution 4: {coords[3]}")
    print("\nThe sum of these coordinates is calculated as follows:")
    # The formatting ensures each number is clearly visible in the final equation.
    print(f"({coords[0]}) + ({coords[1]}) + ({coords[2]}) + ({coords[3]}) = {total_sum}")

solve_and_sum_first_coordinates()