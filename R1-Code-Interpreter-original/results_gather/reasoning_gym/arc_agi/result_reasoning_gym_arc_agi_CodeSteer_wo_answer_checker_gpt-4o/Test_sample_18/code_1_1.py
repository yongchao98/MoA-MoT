def extract_pattern(input_grid):
    # Find the first occurrence of '7' and extract a 3x3 pattern
    for i in range(len(input_grid) - 2):
        for j in range(len(input_grid[i]) - 2):
            if input_grid[i][j] == 7:
                # Extract the 3x3 pattern starting from (i, j)
                pattern = [row[j:j+3] for row in input_grid[i:i+3]]
                return pattern
    return None

def construct_output_grid(pattern, rows, cols):
    # Construct the output grid by repeating the pattern
    output_grid = []
    for i in range(rows):
        for row in pattern:
            output_grid.append(row)
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 2, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 2, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [7, 7, 7, 9, 9, 9, 9, 9, 9],
    [7, 9, 7, 9, 9, 9, 9, 9, 9],
    [9, 7, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9]
]

# Extract the pattern
pattern = extract_pattern(input_grid)

# Construct the output grid
output_grid = construct_output_grid(pattern, 3, 3)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))