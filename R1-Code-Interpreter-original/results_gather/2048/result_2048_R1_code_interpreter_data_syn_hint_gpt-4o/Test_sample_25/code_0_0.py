# Initial grid
grid = [
    [32, 8, 32],
    [0, 32, 16],
    [0, 0, 2]
]

# Initial position of the number to move
pos = (2, 2)

# Sequence of moves
moves = 'uul'

# Function to perform the move
def move_number(grid, pos, moves):
    x, y = pos
    for move in moves:
        if move == 'u':
            if x > 0:
                if grid[x-1][y] == grid[x][y]:
                    grid[x-1][y] *= 2
                    grid[x][y] = 0
                elif grid[x-1][y] == 0:
                    grid[x-1][y] = grid[x][y]
                    grid[x][y] = 0
                x -= 1
        elif move == 'd':
            if x < 2:
                if grid[x+1][y] == grid[x][y]:
                    grid[x+1][y] *= 2
                    grid[x][y] = 0
                elif grid[x+1][y] == 0:
                    grid[x+1][y] = grid[x][y]
                    grid[x][y] = 0
                x += 1
        elif move == 'l':
            if y > 0:
                if grid[x][y-1] == grid[x][y]:
                    grid[x][y-1] *= 2
                    grid[x][y] = 0
                elif grid[x][y-1] == 0:
                    grid[x][y-1] = grid[x][y]
                    grid[x][y] = 0
                y -= 1
        elif move == 'r':
            if y < 2:
                if grid[x][y+1] == grid[x][y]:
                    grid[x][y+1] *= 2
                    grid[x][y] = 0
                elif grid[x][y+1] == 0:
                    grid[x][y+1] = grid[x][y]
                    grid[x][y] = 0
                y += 1
    return grid

# Perform the moves
result_grid = move_number(grid, pos, moves)

# Print the final grid
print(result_grid)