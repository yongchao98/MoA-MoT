# Initial grid
grid = [
    [16, 32],
    [8, 16]
]

# Move sequence
sequence = 'rd'

# Function to move the number in the grid
def move_number(grid, sequence):
    # Start position
    x, y = 0, 0  # (1, 1) in 0-indexed is (0, 0)
    
    for move in sequence:
        if move == 'r':  # Move right
            if y + 1 < 2:  # Check bounds
                y += 1
        elif move == 'd':  # Move down
            if x + 1 < 2:  # Check bounds
                x += 1
        
        # Check for combination
        if grid[x][y] == grid[x-1][y] or grid[x][y] == grid[x][y-1]:
            grid[x][y] *= 2
            if move == 'r':
                grid[x][y-1] = 0
            elif move == 'd':
                grid[x-1][y] = 0
    
    # Set the initial position to 0 after moving
    grid[0][0] = 0
    
    return grid

# Apply the move sequence
result_grid = move_number(grid, sequence)

# Print the result
print(result_grid)