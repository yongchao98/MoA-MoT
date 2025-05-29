# Initial grid
grid = [
    [8, 0, 0, 0, 32],
    [32, 0, 8, 0, 0],
    [0, 0, 0, 32, 2],
    [0, 32, 32, 0, 0],
    [0, 0, 2, 0, 0]
]

# Initial position of the number
x, y = 4, 2  # 0-indexed position (5, 3)

# Move sequence
moves = 'rllddl'

# Function to move the number in the grid
def move_number(grid, x, y, moves):
    for move in moves:
        if move == 'r' and y < 4:
            if grid[x][y + 1] == grid[x][y]:
                grid[x][y + 1] *= 2
                grid[x][y] = 0
            elif grid[x][y + 1] == 0:
                grid[x][y + 1] = grid[x][y]
                grid[x][y] = 0
            y += 1
        elif move == 'l' and y > 0:
            if grid[x][y - 1] == grid[x][y]:
                grid[x][y - 1] *= 2
                grid[x][y] = 0
            elif grid[x][y - 1] == 0:
                grid[x][y - 1] = grid[x][y]
                grid[x][y] = 0
            y -= 1
        elif move == 'd' and x < 4:
            if grid[x + 1][y] == grid[x][y]:
                grid[x + 1][y] *= 2
                grid[x][y] = 0
            elif grid[x + 1][y] == 0:
                grid[x + 1][y] = grid[x][y]
                grid[x][y] = 0
            x += 1
        elif move == 'u' and x > 0:
            if grid[x - 1][y] == grid[x][y]:
                grid[x - 1][y] *= 2
                grid[x][y] = 0
            elif grid[x - 1][y] == 0:
                grid[x - 1][y] = grid[x][y]
                grid[x][y] = 0
            x -= 1
    return grid

# Apply the moves
final_grid = move_number(grid, x, y, moves)

# Print the final grid
print(final_grid)