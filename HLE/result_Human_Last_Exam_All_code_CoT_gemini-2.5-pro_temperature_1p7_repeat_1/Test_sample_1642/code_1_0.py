# Define the size of the chessboard
n = 16

# This problem is known as the Peaceful Queens Problem.
# For a board of size n x n, where n is a multiple of 4,
# the maximum number of peaceful queen pairs (m) is given by the formula:
# m = n * (n - 2) / 4

# Calculate the intermediate term
n_minus_2 = n - 2

# Calculate the numerator
numerator = n * n_minus_2

# Calculate the final value for m
m = numerator // 4

# Print the full equation as requested, showing the intermediate steps
print(f"The calculation for the maximum number m is:")
print(f"{n} * ({n} - 2) / 4 = {n} * {n_minus_2} / 4 = {numerator} / 4 = {m}")
