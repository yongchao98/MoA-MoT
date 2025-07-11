# Step 1: Define the lengths from the initial structure.
initial_stem = 40
initial_branch = 20

# The problem shows a specific path highlighted in white on the 2nd iteration structure.
# We need to calculate the length of this path by summing its constituent parts.

# Step 2: Identify the segments of the white path based on the fractal generation rule.
# The path consists of the main stem, the stem from the first iteration,
# the stem from the second iteration, and the two branches from the second iteration.

# Length of the main stem (from the original structure).
segment1_length = initial_stem

# Length of the stem from the 1st iteration 'Y'. It replaced the initial branch.
segment2_length = initial_branch

# The scaling factor between a stem and its branch is initial_branch / initial_stem = 0.5
scaling_factor = initial_branch / initial_stem

# Length of the stem from the 2nd iteration 'Y'. It replaced a branch from the 1st iteration 'Y'.
# The branch length of the 1st iteration 'Y' was segment2_length * scaling_factor.
segment3_length = segment2_length * scaling_factor

# Length of the two final branches. They are the branches of the 2nd iteration 'Y'.
# Their length is the stem length (segment3_length) times the scaling factor.
segment4_length = segment3_length * scaling_factor
segment5_length = segment3_length * scaling_factor

# Step 3: Calculate the total length by summing the segments.
total_length = segment1_length + segment2_length + segment3_length + segment4_length + segment5_length

# Step 4: Print the calculation step-by-step.
print("The total length of the white path is calculated by summing its segments:")
print(f"Segment 1 (Initial Stem): {segment1_length}")
print(f"Segment 2 (1st Iteration Stem): {segment2_length}")
print(f"Segment 3 (2nd Iteration Stem): {segment3_length}")
print(f"Segments 4 & 5 (2nd Iteration Branches): {segment4_length} and {segment5_length}")
print("\nFinal Equation:")
print(f"{int(segment1_length)} + {int(segment2_length)} + {int(segment3_length)} + {int(segment4_length)} + {int(segment5_length)} = {int(total_length)}")
