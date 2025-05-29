# Initial grid
grid = [
    [4, 8, 0, 8],
    [16, 0, 0, 0],
    [16, 2, 0, 0],
    [0, 32, 16, 0]
]

# Initial position of the number to move
position = (3, 0)  # 0-indexed, corresponds to (4, 1) in 1-indexed

# Sequence of moves
moves = 'udul'

# Function to move the number
def move_number(grid, position, moves):
    for move in moves:
        x, y = position
        if move == 'u' and x > 0:
            position = (x - 1, y)
        elif move == 'd' and x < len(grid) - 1:
            position = (x + 1, y)
        elif move == 'l' and y > 0:
            position = (x, y - 1)
        elif move == 'r' and y < len(grid[0]) - 1:
            position = (x, y + 1)
        
        # Check for combination
        if grid[position[0]][position[1]] == grid[x][y]:
            grid[position[0]][position[1]] *= 2
            grid[x][y] = 0

    return grid

# Apply the moves
result_grid = move_number(grid, position, moves)

# Print the resulting grid
print(result_grid)