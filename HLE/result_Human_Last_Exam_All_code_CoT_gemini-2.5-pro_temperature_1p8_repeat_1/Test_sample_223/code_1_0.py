# The total number of chairs in the row.
N = 20

# The formula for the maximum number of occupants in this puzzle is
# known to be floor(2 * N / 3). This is derived from analyzing the
# stable configurations where no more people can be added to increase the count.
# The most efficient stable patterns are arrangements of 'OOX' blocks.

# Constants for the formula
multiplier = 2
divisor = 3

# Calculate the maximum number of occupied chairs using integer division for floor.
max_occupied = (multiplier * N) // divisor

# As requested, the output shows the calculation with all numbers involved.
print(f"For {N} chairs, the maximum number of occupied chairs can be found using the formula: floor(2 * N / 3)")
print(f"Plugging in the numbers:")
print(f"Maximum = floor({multiplier} * {N} / {divisor}) = {max_occupied}")