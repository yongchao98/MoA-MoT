# Initial grid
grid = [
    [0, 0, 0, 16, 0],
    [0, 0, 0, 16, 0],
    [8, 0, 2, 0, 0],
    [0, 0, 0, 0, 0],
    [2, 0, 8, 32, 0]
]

# Initial position of the number 8
position = (4, 2)  # (row, column) in 0-indexed

# Movement sequence
moves = 'llruu'

# Function to move the number according to the sequence
def move_number(grid, position, moves):
    for move in moves:
        row, col = position
        if move == 'l':
            new_col = col - 1
            if new_col >= 0 and (grid[row][new_col] == 0 or grid[row][new_col] == grid[row][col]):
                grid[row][new_col] += grid[row][col]
                grid[row][col] = 0
                position = (row, new_col)
        elif move == 'r':
            new_col = col + 1
            if new_col < 5 and (grid[row][new_col] == 0 or grid[row][new_col] == grid[row][col]):
                grid[row][new_col] += grid[row][col]
                grid[row][col] = 0
                position = (row, new_col)
        elif move == 'u':
            new_row = row - 1
            if new_row >= 0 and (grid[new_row][col] == 0 or grid[new_row][col] == grid[row][col]):
                grid[new_row][col] += grid[row][col]
                grid[row][col] = 0
                position = (new_row, col)
        elif move == 'd':
            new_row = row + 1
            if new_row < 5 and (grid[new_row][col] == 0 or grid[new_row][col] == grid[row][col]):
                grid[new_row][col] += grid[row][col]
                grid[row][col] = 0
                position = (new_row, col)
    return grid

# Apply the moves
final_grid = move_number(grid, position, moves)

# Print the final grid
print(final_grid)