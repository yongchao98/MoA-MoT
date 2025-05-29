import numpy as np

# Initial counts of each block type
initial_blocks = np.array([5, 4, 4, 0, 0, 0, 0, 0, 0])

# Transformation matrix representing the rules
# Each row corresponds to a rule, and each column corresponds to a block type
transformation_matrix = np.array([
    [-1, -1, -1, 1, 0, 0, 0, 0, 0],  # Rule 1
    [-1, -1, 0, 0, 0, 1, 0, 0, 0],   # Rule 2
    [0, -1, -1, 0, 1, 0, 0, 0, 0],   # Rule 3
    [0, 0, -2, 0, 0, 1, 0, 0, 0],    # Rule 4
    [0, 0, 0, -1, 0, -1, 1, 1, 0],   # Rule 5
    [0, 0, 0, 0, -2, 0, 0, 0, 1]     # Rule 6
])

# Function to apply transformations
def apply_transformations(blocks, transformation_matrix):
    changes = True
    while changes:
        changes = False
        for i, transformation in enumerate(transformation_matrix):
            # Calculate the maximum number of times the transformation can be applied
            max_applications = min(
                [blocks[j] // -transformation[j] if transformation[j] < 0 else float('inf') for j in range(len(blocks))]
            )
            if max_applications > 0:
                # Apply the transformation
                blocks += max_applications * transformation
                changes = True
    return blocks

# Apply the transformations
final_blocks = apply_transformations(initial_blocks, transformation_matrix)

# Prepare the final answer
final_answer = ''.join(map(str, final_blocks))

# Output the final answer
print(f"Your answer: {final_answer}")