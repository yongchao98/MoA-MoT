# Initial grid
grid = [
    [0, 0, 0, 32, 0],
    [0, 0, 0, 0, 16],
    [0, 2, 16, 0, 32],
    [0, 0, 2, 2, 0],
    [0, 32, 16, 0, 0]
]

# Initial position of the number to move
x, y = 2, 4  # (3, 5) in 0-indexed

# Movement sequence
moves = 'rulrdr'

# Function to move the number in the grid
def move_number(grid, x, y, moves):
    for move in moves:
        if move == 'r':  # Move right
            if y < 4:
                if grid[x][y+1] == grid[x][y]:
                    grid[x][y+1] *= 2
                    grid[x][y] = 0
                elif grid[x][y+1] == 0:
                    grid[x][y+1] = grid[x][y]
                    grid[x][y] = 0
                y = min(y + 1, 4)
        elif move == 'u':  # Move up
            if x > 0:
                if grid[x-1][y] == grid[x][y]:
                    grid[x-1][y] *= 2
                    grid[x][y] = 0
                elif grid[x-1][y] == 0:
                    grid[x-1][y] = grid[x][y]
                    grid[x][y] = 0
                x = max(x - 1, 0)
        elif move == 'l':  # Move left
            if y > 0:
                if grid[x][y-1] == grid[x][y]:
                    grid[x][y-1] *= 2
                    grid[x][y] = 0
                elif grid[x][y-1] == 0:
                    grid[x][y-1] = grid[x][y]
                    grid[x][y] = 0
                y = max(y - 1, 0)
        elif move == 'd':  # Move down
            if x < 4:
                if grid[x+1][y] == grid[x][y]:
                    grid[x+1][y] *= 2
                    grid[x][y] = 0
                elif grid[x+1][y] == 0:
                    grid[x+1][y] = grid[x][y]
                    grid[x][y] = 0
                x = min(x + 1, 4)
    return grid

# Apply the moves
result_grid = move_number(grid, x, y, moves)

# Print the resulting grid
print(result_grid)