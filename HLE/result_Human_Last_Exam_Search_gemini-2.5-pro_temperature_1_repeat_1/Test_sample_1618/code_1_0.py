import numpy as np

# Define the matrix P
P = np.array([
    [0.985, 0.015, 0, 0],
    [0.5, 0.4, 0.1, 0],
    [0, 0.99, 0, 0.1],
    [0, 1, 0, 0]
])

# Define the power n
n = 3431

# Calculate P raised to the power of n
P_n = np.linalg.matrix_power(P, n)

# Calculate the sum of the squares of the elements of P_n
sum_of_squares = np.sum(np.square(P_n))

# The result is then rounded to three decimal places.