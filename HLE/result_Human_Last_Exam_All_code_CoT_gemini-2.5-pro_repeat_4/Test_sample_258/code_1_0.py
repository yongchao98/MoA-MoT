# The radius of the circle.
R = 500

# The number of cells a circle crosses depends on the position of its center.
# For a large integer radius R, known results from number theory provide the
# formulas for the minimum and maximum number of crossed cells.

# The maximum number of cells is achieved when the center of the circle is
# placed at (k + 0.5, l + 0.5) for any integers k, l.
# The formula is N_max = 8 * R - 4.
max_cells = 8 * R - 4

# The minimum number of cells is achieved when the center of the circle is
# placed at (k + epsilon, l + epsilon) for any integers k, l and a very
# small epsilon > 0.
# The formula is N_min = 8 * R - 8.
min_cells = 8 * R - 8

# Print the calculations and the final result.
print(f"The radius of the circle is R = {R}.")
print(f"The maximal number of cells is calculated as 8 * {R} - 4 = {max_cells}.")
print(f"The minimal number of cells is calculated as 8 * {R} - 8 = {min_cells}.")
print(f"The answer in the form (min, max) is ({min_cells}, {max_cells}).")
