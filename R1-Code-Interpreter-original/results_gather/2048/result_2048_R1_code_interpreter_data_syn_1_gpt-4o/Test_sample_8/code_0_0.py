# Initial grid
grid = [
    [0, 4],
    [0, 4]
]

# Move sequence
moves = 'dr'

# Function to perform the moves
def move_2048(grid, moves):
    for move in moves:
        if move == 'd':  # Move down
            if grid[0][1] != 0 and grid[1][1] == 0:
                grid[1][1] = grid[0][1]
                grid[0][1] = 0
            elif grid[0][1] == grid[1][1]:
                grid[1][1] *= 2
                grid[0][1] = 0
        elif move == 'r':  # Move right
            # Since the number is already at the rightmost column, no change
            pass
    return grid

# Perform the moves
result_grid = move_2048(grid, moves)

# Print the result
print(result_grid)