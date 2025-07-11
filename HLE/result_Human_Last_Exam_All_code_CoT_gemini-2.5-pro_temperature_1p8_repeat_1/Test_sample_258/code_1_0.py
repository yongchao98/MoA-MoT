# The radius of the circle
R = 500

# The number of grid line crossings can be calculated as follows.
# Number of vertical lines crossed (K_x):
# The lines x=n are crossed if x_0 - R < n < x_0 + R.
# Since R is an integer and x_0 is not, the number of integers in this interval of length 2R=1000 is always 1000.
K_x = 2 * R
# Number of horizontal lines crossed (K_y):
K_y = 2 * R

# Total number of times the circle crosses grid lines:
# Each vertical line is crossed twice, and each horizontal line is crossed twice.
num_crossings = 2 * K_x + 2 * K_y

# Based on geometric considerations, the number of cells crossed (C) can vary slightly.
# The maximum number of cells crossed is equal to the number of line crossings.
# This occurs when the circle's center is infinitesimally close to a grid intersection (e.g., (epsilon, epsilon)).
max_cells = 8 * R

# The minimum number of cells crossed occurs when the center is at (I+0.5, J+0.5).
# In this highly symmetric case, 4 cells are effectively not counted, leading to a reduction.
min_cells = 8 * R - 4

# Print the final result
print(f"R = {R}")
print(f"Minimal number of cells crossed = 8 * R - 4 = 8 * {R} - 4 = {min_cells}")
print(f"Maximal number of cells crossed = 8 * R = 8 * {R} = {max_cells}")
print(f"The answer is ({min_cells}, {max_cells})")
