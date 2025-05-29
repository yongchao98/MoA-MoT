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
                if grid[x][y] == grid[x][y + 1]:  # Check for combination
                    grid[x][y + 1] *= 2
                    grid[x][y] = 0
                else:
                    grid[x][y + 1] = grid[x][y]
                    grid[x][y] = 0
                y += 1
        elif move == 'd':  # Move down
            if x + 1 < 2:  # Check bounds
                if grid[x][y] == grid[x + 1][y]:  # Check for combination
                    grid[x + 1][y] *= 2
                    grid[x][y] = 0
                else:
                    grid[x + 1][y] = grid[x][y]
                    grid[x][y] = 0
                x += 1
    
    return grid

# Apply the move sequence
result_grid = move_number(grid, sequence)

# Print the result
print(result_grid)