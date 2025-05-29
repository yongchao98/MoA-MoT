def transform_grid(input_grid):
    # Flatten the grid to count occurrences of each number
    from collections import Counter
    flat_grid = [num for row in input_grid for num in row]
    count = Counter(flat_grid)
    
    # Find the majority number
    majority_number = count.most_common(1)[0][0]
    
    # Transform the grid
    output_grid = [[majority_number for _ in row] for row in input_grid]
    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6],
    [6, 6, 6],
    [6, 6, 6],
    [6, 6, 6],
    [6, 6, 6]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))