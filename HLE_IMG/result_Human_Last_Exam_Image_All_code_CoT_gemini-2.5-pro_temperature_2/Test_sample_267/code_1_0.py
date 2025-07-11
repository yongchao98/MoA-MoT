# Step 1: Define the lengths from the initial structure (Iteration 0).
main_trunk = 40
initial_branch = 20

# Step 2: Determine the scaling factor from the initial structure.
# The ratio of a new branch to its trunk is constant.
scaling_factor = initial_branch / main_trunk

# Step 3: Identify the segments of the white path based on the fractal generation rule.
# The path follows the main trunk and then the right-most branch at each iteration.

# Length of the main trunk (the first white segment).
l1 = main_trunk

# Length of the right branch from Iteration 0 (the second white segment).
l2 = initial_branch

# Length of the right branch from Iteration 1 (the third white segment).
# This is calculated by scaling the previous branch length.
l3 = l2 * scaling_factor

# Lengths of the two final branches from Iteration 2 (the fourth and fifth white segments).
# These are calculated by scaling the previous branch length.
l4 = l3 * scaling_factor
l5 = l3 * scaling_factor

# Step 4: Calculate the total length of the white branches.
total_length = l1 + l2 + l3 + l4 + l5

# Print the breakdown of the calculation and the final answer.
# Note: We use int() to display the lengths as whole numbers.
print(f"The lengths of the segments that form the white path are {int(l1)}, {int(l2)}, {int(l3)}, {int(l4)}, and {int(l5)}.")
print(f"The total length is the sum of these segments:")
print(f"{int(l1)} + {int(l2)} + {int(l3)} + {int(l4)} + {int(l5)} = {int(total_length)}")