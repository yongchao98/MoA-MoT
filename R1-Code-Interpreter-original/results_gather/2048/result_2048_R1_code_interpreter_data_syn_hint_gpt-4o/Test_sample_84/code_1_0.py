# Initial grid
grid = [
    [32, 0, 32, 8, 0],
    [16, 0, 0, 0, 0],
    [0, 0, 16, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 2, 0, 0, 0]
]

# Initial position of the number
position = (4, 1)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'llldud'

# Function to move the number
def move_number(grid, position, moves):
    for move in moves:
        row, col = position
        if move == 'l':
            if col > 0:
                if grid[row][col - 1] == grid[row][col]:
                    grid[row][col - 1] *= 2
                    grid[row][col] = 0
                elif grid[row][col - 1] == 0:
                    grid[row][col - 1] = grid[row][col]
                    grid[row][col] = 0
                position = (row, col - 1)
        elif move == 'r':
            if col < 4:
                if grid[row][col + 1] == grid[row][col]:
                    grid[row][col + 1] *= 2
                    grid[row][col] = 0
                elif grid[row][col + 1] == 0:
                    grid[row][col + 1] = grid[row][col]
                    grid[row][col] = 0
                position = (row, col + 1)
        elif move == 'u':
            if row > 0:
                if grid[row - 1][col] == grid[row][col]:
                    grid[row - 1][col] *= 2
                    grid[row][col] = 0
                elif grid[row - 1][col] == 0:
                    grid[row - 1][col] = grid[row][col]
                    grid[row][col] = 0
                position = (row - 1, col)
        elif move == 'd':
            if row < 4:
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