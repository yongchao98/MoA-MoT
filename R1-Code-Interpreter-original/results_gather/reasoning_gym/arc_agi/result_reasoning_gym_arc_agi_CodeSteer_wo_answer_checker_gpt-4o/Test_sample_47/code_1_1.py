def transform_grid(input_grid):
    # Define the regions in the input grid that correspond to each cell in the output grid
    regions = [
        [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)],  # Top-left
        [(0, 4), (0, 5), (0, 6), (1, 4), (1, 5), (1, 6), (2, 4), (2, 5), (2, 6)],  # Top-middle
        [(0, 8), (0, 9), (0, 10), (1, 8), (1, 9), (1, 10), (2, 8), (2, 9), (2, 10)],  # Top-right
        [(4, 0), (4, 1), (4, 2), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (6, 2)],  # Middle-left
        [(4, 4), (4, 5), (4, 6), (5, 4), (5, 5), (5, 6), (6, 4), (6, 5), (6, 6)],  # Center
        [(4, 8), (4, 9), (4, 10), (5, 8), (5, 9), (5, 10), (6, 8), (6, 9), (6, 10)],  # Middle-right
        [(8, 0), (8, 1), (8, 2), (9, 0), (9, 1), (9, 2), (10, 0), (10, 1), (10, 2)],  # Bottom-left
        [(8, 4), (8, 5), (8, 6), (9, 4), (9, 5), (9, 6), (10, 4), (10, 5), (10, 6)],  # Bottom-middle
        [(8, 8), (8, 9), (8, 10), (9, 8), (9, 9), (9, 10), (10, 8), (10, 9), (10, 10)]  # Bottom-right
    ]
    
    # Initialize the output grid
    output_grid = []
    
    # Process each region
    for region in regions:
        # Count the number of '1's in the region
        count_ones = sum(1 for (r, c) in region if input_grid[r][c] == 1)
        
        # Determine the output value for this region
        output_value = 4 if count_ones > 0 else 6
        output_grid.append(output_value)
    
    # Reshape the output grid into a 3x3 matrix
    output_grid = [output_grid[i:i+3] for i in range(0, 9, 3)]
    
    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 9, 6, 6, 6, 9, 6, 1, 6],
    [6, 1, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [6, 6, 6, 9, 6, 1, 6, 9, 1, 6, 6],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 9, 6, 6, 1, 9, 6, 6, 6],
    [6, 6, 6, 9, 6, 1, 6, 9, 6, 6, 1],
    [6, 1, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [1, 6, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [6, 6, 6, 9, 1, 6, 6, 9, 6, 1, 6],
    [1, 6, 6, 9, 6, 6, 1, 9, 6, 6, 6]
]

output_grid = transform_grid(input_grid)
print("<<<", output_grid, ">>>")