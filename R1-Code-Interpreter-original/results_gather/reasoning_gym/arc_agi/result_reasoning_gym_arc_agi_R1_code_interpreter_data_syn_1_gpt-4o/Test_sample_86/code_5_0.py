def transform_grid(input_grid):
    # Initialize the output grid with 8s
    output_grid = [[8] * 6 for _ in range(3)]
    
    # Function to find sequences of 9s in a row
    def find_nine_sequences(grid):
        sequences = []
        for row in grid:
            start = None
            for i, val in enumerate(row):
                if val == 9 and start is None:
                    start = i
                elif val != 9 and start is not None:
                    sequences.append((start, i - 1))
                    start = None
            if start is not None:
                sequences.append((start, len(row) - 1))
        return sequences
    
    # Find sequences of 9s in the entire grid
    sequences = find_nine_sequences(input_grid)
    
    if sequences:
        first_start, _ = sequences[0]
        _, last_end = sequences[-1]
        if first_start < 6:
            output_grid[0][first_start] = 1
            output_grid[2][first_start] = 1
        if last_end < 6:
            output_grid[0][last_end] = 3
            output_grid[2][last_end] = 3
    
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 8, 8, 8, 8, 8, 8],
    [8, 9, 8, 8, 8, 8, 8, 9, 8],
    [8, 8, 8, 9, 8, 8, 9, 9, 9],
    [8, 8, 9, 9, 9, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 9, 8, 8, 8],
    [9, 9, 9, 8, 9, 9, 9, 8, 8],
    [8, 9, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 9, 9, 9],
    [8, 8, 9, 9, 9, 8, 8, 9, 8],
    [8, 8, 8, 9, 8, 8, 8, 8, 8]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))