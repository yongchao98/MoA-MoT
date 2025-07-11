import math

# Step 1: Define the number of simple paths.
# A path can go from A to P, then take one of the three routes from P to Q,
# and finally go from Q to B. These are the simple, non-self-intersecting paths.
num_simple_paths = 3

# Step 2: Define the number of looping paths.
# A path is allowed to self-intersect. This means it can loop.
# A path can go from P to Q (choice 1), then back from Q to P (choice 2),
# and finally from P to Q again (choice 3).
# This creates a path that uses all three P-Q edges in a specific order.
# The number of ways to order these three distinct edges is 3 factorial.
num_looping_paths = math.factorial(3)

# Step 3: Calculate the total number of paths.
# The total number of distinct paths (trails) is the sum of the simple paths
# and the looping paths.
total_paths = num_simple_paths + num_looping_paths

# Step 4: Print the final equation and the result.
# The problem asks for the final equation.
print(f"The number of paths can be calculated by summing the simple paths and the looping paths.")
print(f"Number of simple paths: {num_simple_paths}")
print(f"Number of looping paths (permutations of 3 edges): {num_looping_paths}")
print(f"Total paths = {num_simple_paths} + {num_looping_paths} = {total_paths}")