# The radius of the circle
R = 500

# Step 1: Calculate the number of horizontal and vertical grid lines crossed.
# The number of horizontal lines crossed (N_h) is the number of integers in the
# interval (yc - R, yc + R). Since the length of this interval is 2*R = 1000 and
# its endpoints are not integers (due to the non-tangency condition),
# N_h is always 1000.
N_h = 2 * R
# Similarly, the number of vertical lines crossed (N_v) is 1000.
N_v = 2 * R

# Step 2: Calculate the total number of grid line crossings.
# Each line is crossed twice by the circle.
num_crossings = 2 * N_h + 2 * N_v

# Step 3: Determine the minimal and maximal number of cells.
# The number of cells crossed can be at most the number of crossings.
# This maximum is achieved when the center of the circle is placed near
# (integer + 0.5, integer + 0.5), which ensures no cell re-entries.
max_cells = num_crossings

# The number of cells can be reduced by up to 4 due to the "poke" effect
# at the four cardinal extrema of the circle (North, South, East, West).
# This minimum is achieved when the center is placed very close to a
# grid intersection (e.g., (integer + epsilon, integer + epsilon)).
# This creates 4 cell re-entries, reducing the count by 4.
min_cells = num_crossings - 4

# The final answer is the tuple (min_cells, max_cells).
# We print the calculation for clarity.
print(f"Radius (R): {R}")
print(f"Number of horizontal lines crossed (N_h): {N_h}")
print(f"Number of vertical lines crossed (N_v): {N_v}")
print(f"Total grid line crossings: 2 * (N_h + N_v) = 2 * ({N_h} + {N_v}) = {num_crossings}")
print(f"Maximal number of cells = Total crossings = {max_cells}")
print(f"Minimal number of cells = Total crossings - 4 = {num_crossings} - 4 = {min_cells}")
print(f"The minimal and maximal numbers of cells are: ({min_cells}, {max_cells})")
