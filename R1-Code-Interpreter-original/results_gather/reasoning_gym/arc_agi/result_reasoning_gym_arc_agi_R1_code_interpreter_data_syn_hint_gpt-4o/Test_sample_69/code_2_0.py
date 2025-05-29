def extract_output_grid(input_grid):
    # Convert the input grid into a 2D list
    grid = [list(map(int, input_grid[i:i+21])) for i in range(0, len(input_grid), 21)]
    
    # Identify key numbers based on the examples
    key_numbers = [5, 4, 6, 2, 1]  # These numbers frequently appear in the output grids
    
    # Initialize the output grid
    output_grid = []
    
    # Extract rows and columns based on the presence of key numbers
    for row in grid:
        if any(num in row for num in key_numbers):
            # Extract columns where key numbers appear
            selected_row = [row[i] for i in range(len(row)) if row[i] in key_numbers]
            if selected_row:
                output_grid.append(selected_row)
    
    # Limit the output grid to a specific size if necessary
    # Here, we assume a 3x3 grid based on the examples
    output_grid = [row[:3] for row in output_grid[:3]]
    
    return output_grid

# Test input grid
input_grid = [
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3,
    3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3,
    3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3,
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3,
    3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3,
    3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3,
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3,
    3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3,
    3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3,
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3,
    3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3,
    3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3,
    3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3
]

output_grid = extract_output_grid(input_grid)
print(output_grid)