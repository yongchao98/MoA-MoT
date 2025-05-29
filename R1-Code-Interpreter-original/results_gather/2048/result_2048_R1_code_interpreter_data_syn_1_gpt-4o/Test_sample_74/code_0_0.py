# Initial grid
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 2, 0, 32],
    [0, 0, 0, 0, 0],
    [0, 0, 16, 0, 0],
    [16, 2, 0, 0, 0]
]

# Initial position of the number
position = (1, 4)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'rrurl'

# Function to move the number in the grid
def move_number(grid, position, moves):
    for move in moves:
        row, col = position
        if move == 'r' and col < 4:
            # Move right
            if grid[row][col + 1] == grid[row][col]:
                grid[row][col + 1] *= 2
                grid[row][col] = 0
            elif grid[row][col + 1] == 0:
                grid[row][col + 1] = grid[row][col]
                grid[row][col] = 0
            position = (row, col + 1)
        elif move == 'l' and col > 0:
            # Move left
            if grid[row][col - 1] == grid[row][col]:
                grid[row][col - 1] *= 2
                grid[row][col] = 0
            elif grid[row][col - 1] == 0:
                grid[row][col - 1] = grid[row][col]
                grid[row][col] = 0
            position = (row, col - 1)
        elif move == 'u' and row > 0:
            # Move up
            if grid[row - 1][col] == grid[row][col]:
                grid[row - 1][col] *= 2
                grid[row][col] = 0
            elif grid[row - 1][col] == 0:
                grid[row - 1][col] = grid[row][col]
                grid[row][col] = 0
            position = (row - 1, col)
        elif move == 'd' and row < 4:
            # Move down
            if grid[row + 1][col] == grid[row][col]:
                grid[row + 1][col] *= 2
                grid[row][col] = 0
            elif grid[row + 1][col] == 0:
                grid[row + 1][col] = grid[row][col]
                grid[row][col] = 0
            position = (row + 1, col)
    return grid

# Apply the moves
result_grid = move_number(grid, position, moves)

# Print the resulting grid
print(result_grid)