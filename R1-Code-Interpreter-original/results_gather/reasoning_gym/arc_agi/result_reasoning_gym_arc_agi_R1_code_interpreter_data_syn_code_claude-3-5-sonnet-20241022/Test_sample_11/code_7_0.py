def process_grid(input_grid):
    # Create a deep copy of the input grid using list comprehension
    output_grid = [[num for num in row] for row in input_grid]
    
    # Get dimensions
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Process each cell
    for i in range(rows):
        for j in range(cols):
            # Check if current cell contains 0
            if input_grid[i][j] == 0:
                # Explicitly check row position and replace accordingly
                if i <= 4:  # Upper half (rows 0-4)
                    output_grid[i][j] = input_grid[0][j]  # Use number from first row
                else:  # Lower half (rows 5-9)
                    output_grid[i][j] = input_grid[rows-1][j]  # Use number from last row
    
    # Convert to string format with proper spacing
    result = ""
    for row in output_grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())  # For verification
    return output_grid

# Test input
test_grid = [
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,0,8,8,8,0,8,8,6],
    [4,8,8,8,0,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,0,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,8,8,0,8,8,0,6],
    [4,0,8,8,8,8,8,8,8,6]
]

process_grid(test_grid)