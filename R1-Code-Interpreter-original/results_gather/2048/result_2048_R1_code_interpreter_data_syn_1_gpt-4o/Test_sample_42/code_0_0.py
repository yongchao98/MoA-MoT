def move_up(grid):
    for col in range(4):
        for row in range(1, 4):
            if grid[row][col] != 0:
                for move_row in range(row, 0, -1):
                    if grid[move_row-1][col] == 0:
                        grid[move_row-1][col] = grid[move_row][col]
                        grid[move_row][col] = 0
                    elif grid[move_row-1][col] == grid[move_row][col]:
                        grid[move_row-1][col] *= 2
                        grid[move_row][col] = 0
                        break
                    else:
                        break
    return grid

def move_down(grid):
    for col in range(4):
        for row in range(2, -1, -1):
            if grid[row][col] != 0:
                for move_row in range(row, 3):
                    if grid[move_row+1][col] == 0:
                        grid[move_row+1][col] = grid[move_row][col]
                        grid[move_row][col] = 0
                    elif grid[move_row+1][col] == grid[move_row][col]:
                        grid[move_row+1][col] *= 2
                        grid[move_row][col] = 0
                        break
                    else:
                        break
    return grid

def move_left(grid):
    for row in range(4):
        for col in range(1, 4):
            if grid[row][col] != 0:
                for move_col in range(col, 0, -1):
                    if grid[row][move_col-1] == 0:
                        grid[row][move_col-1] = grid[row][move_col]
                        grid[row][move_col] = 0
                    elif grid[row][move_col-1] == grid[row][move_col]:
                        grid[row][move_col-1] *= 2
                        grid[row][move_col] = 0
                        break
                    else:
                        break
    return grid

def move_right(grid):
    for row in range(4):
        for col in range(2, -1, -1):
            if grid[row][col] != 0:
                for move_col in range(col, 3):
                    if grid[row][move_col+1] == 0:
                        grid[row][move_col+1] = grid[row][move_col]
                        grid[row][move_col] = 0
                    elif grid[row][move_col+1] == grid[row][move_col]:
                        grid[row][move_col+1] *= 2
                        grid[row][move_col] = 0
                        break
                    else:
                        break
    return grid

# Initial grid
grid = [
    [0, 0, 4, 0],
    [0, 0, 0, 0],
    [0, 0, 4, 0],
    [2, 0, 16, 16]
]

# Sequence of moves
moves = 'urdu'

# Apply each move
for move in moves:
    if move == 'u':
        grid = move_up(grid)
    elif move == 'd':
        grid = move_down(grid)
    elif move == 'l':
        grid = move_left(grid)
    elif move == 'r':
        grid = move_right(grid)

print(grid)