def transform_grid(input_grid):
    # Convert the input list to a 5x5 grid
    grid = [input_grid[i:i+5] for i in range(0, 25, 5)]
    
    # Initialize the output grid
    output_grid = []
    
    # First row is copied directly
    output_grid.append(grid[0])
    
    # Second row is the first row shifted right by one
    output_grid.append([0] + grid[0][:4])
    
    # Third row is a copy of the first row
    output_grid.append(grid[0])
    
    # Fourth row is the first row shifted right by one
    output_grid.append([0] + grid[0][:4])
    
    # Fifth row is a copy of the first row
    output_grid.append(grid[0])
    
    # Flatten the output grid
    return [num for row in output_grid for num in row]

# Test input
input_grid = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")