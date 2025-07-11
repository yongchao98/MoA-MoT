import math

# Define the parameters from the problem
num_pirates = 9
num_gunpoint_pairs = 16

# The optimal structure is the Wheel Graph W_9.
# It has a central hub and a rim of (num_pirates - 1) vertices.
num_rim_vertices = num_pirates - 1

# 1. Count the cycles on the rim. There is only one.
cycles_on_rim = 1

# 2. Count the cycles that include the central hub.
# This is equivalent to choosing any 2 vertices from the rim.
# We use the combination formula C(n, k) = n! / (k! * (n-k)!).
cycles_with_hub = math.comb(num_rim_vertices, 2)

# The total number of standoffs is the sum of these two counts.
total_standoffs = cycles_on_rim + cycles_with_hub

# Print the final equation and the result.
print(f"The maximum number of standoffs is found by counting the cycles in a Wheel Graph with {num_pirates} vertices.")
print(f"This is the sum of the single outer rim cycle and the cycles formed by choosing 2 out of the {num_rim_vertices} rim vertices.")
print(f"Maximum standoffs = {cycles_on_rim} + C({num_rim_vertices}, 2)")
print(f"Maximum standoffs = {cycles_on_rim} + {cycles_with_hub}")
print(f"Total = {total_standoffs}")
