from collections import Counter

def transform_grid(input_grid):
    # Flatten the grid to count frequencies
    flat_grid = [num for row in input_grid for num in row]
    # Find the most frequent number
    most_common_num, _ = Counter(flat_grid).most_common(1)[0]
    
    # Transform the grid
    output_grid = [[5 if num != most_common_num else num for num in row] for row in input_grid]
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 3, 8],
    [8, 8, 8, 8, 8, 8],
    [8, 8, 3, 8, 3, 8],
    [8, 2, 8, 8, 8, 3],
    [8, 8, 8, 2, 8, 8],
    [8, 8, 8, 8, 8, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))