# Step 1: Define the lengths of the initial structure (Iteration 0)
l0 = 40.0
b0 = 20.0

# Step 2: Calculate the lengths for the structure added in Iteration 1
# This new structure grows on a branch of length b0.
# The new stem (l1) is half the parent branch's length.
# The new branches (b1) are half the new stem's length.
l1 = b0 / 2
b1 = l1 / 2

# Step 3: Calculate the lengths for the structure added in Iteration 2
# This structure grows on a branch from Iteration 1, which has length b1.
l2 = b1 / 2
b2 = l2 / 2

# Step 4: Identify the segments of the white path and sum their lengths
# The white path consists of:
# - The initial stem (l0)
# - The first right branch (b0)
# - The stem from iteration 1 (l1)
# - The left branch from iteration 1 (b1)
# - The entire structure from iteration 2 (stem l2 + two branches b2)

# Calculate the total length of the small Y-structure from iteration 2
path_y2_length = l2 + b2 + b2

# Calculate the total length of the white path
total_length = l0 + b0 + l1 + b1 + path_y2_length

# Step 5: Print the final equation with all its components
print(f"The total length is the sum of the segments: {l0} + {b0} + {l1} + {b1} + ({l2} + {b2} + {b2}) = {total_length}")
