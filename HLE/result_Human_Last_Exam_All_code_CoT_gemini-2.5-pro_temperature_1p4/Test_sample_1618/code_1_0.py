import numpy as np

# Define the matrix P and the power n
P = np.array([
    [0.985, 0.015, 0, 0],
    [0.5, 0.4, 0.1, 0],
    [0, 0.99, 0, 0.1],
    [0, 1, 0, 0]
])
power = 3431

# Calculate P raised to the power of 3431
P_pow = np.linalg.matrix_power(P, power)

# When P is raised to a large power, the resulting matrix has four identical rows.
# We can take the first row as the representative row.
representative_row = P_pow[0]

# Calculate the squares of the elements in the representative row
row_element_squares = representative_row**2

# The total sum of squares is 4 times the sum of squares of one row.
sum_of_squares_row = np.sum(row_element_squares)
total_sum_of_squares = 4 * sum_of_squares_row

# Print the detailed calculation
print(f"The matrix P raised to the power of {power} results in a matrix where each row is identical.")
print(f"The representative row is: [{representative_row[0]:.6f}, {representative_row[1]:.6f}, {representative_row[2]:.6f}, {representative_row[3]:.6f}]")
print("\nThe sum of the squares of the elements is calculated as follows:")
print("4 * (r1^2 + r2^2 + r3^2 + r4^2)")
print(f"= 4 * ({representative_row[0]:.6f}^2 + {representative_row[1]:.6f}^2 + {representative_row[2]:.6f}^2 + {representative_row[3]:.6f}^2)")
print(f"= 4 * ({row_element_squares[0]:.6f} + {row_element_squares[1]:.6f} + {row_element_squares[2]:.6f} + {row_element_squares[3]:.6f})")
print(f"= 4 * ({sum_of_squares_row:.6f})")
print(f"= {total_sum_of_squares:.6f}")
print("\nFinal Answer:")
print(f"The sum of the squares of the elements of P^{power} is {total_sum_of_squares:.3f}")