def move_left(grid):
    for i in range(4):
        new_row = [x for x in grid[i] if x != 0]
        for j in range(1, len(new_row)):
            if new_row[j] == new_row[j - 1]:
                new_row[j - 1] *= 2
                new_row[j] = 0
        new_row = [x for x in new_row if x != 0]
        grid[i] = new_row + [0] * (4 - len(new_row))
    return grid

def move_up(grid):
    for j in range(4):
        new_col = [grid[i][j] for i in range(4) if grid[i][j] != 0]
        for i in range(1, len(new_col)):
            if new_col[i] == new_col[i - 1]:
                new_col[i - 1] *= 2
                new_col[i] = 0
        new_col = [x for x in new_col if x != 0]
        for i in range(4):
            grid[i][j] = new_col[i] if i < len(new_col) else 0
    return grid

def move_down(grid):
    for j in range(4):
        new_col = [grid[i][j] for i in range(4) if grid[i][j] != 0]
        for i in range(len(new_col) - 2, -1, -1):
            if new_col[i] == new_col[i + 1]:
                new_col[i + 1] *= 2
                new_col[i] = 0
        new_col = [x for x in new_col if x != 0]
        for i in range(4):
            grid[3 - i][j] = new_col[len(new_col) - 1 - i] if i < len(new_col) else 0
    return grid

def move_right(grid):
    for i in range(4):
        new_row = [x for x in grid[i] if x != 0]
        for j in range(len(new_row) - 2, -1, -1):
            if new_row[j] == new_row[j + 1]:
                new_row[j + 1] *= 2
                new_row[j] = 0
        new_row = [x for x in new_row if x != 0]
        grid[i] = [0] * (4 - len(new_row)) + new_row
    return grid

def apply_moves(grid, moves):
    for move in moves:
        if move == 'l':
            grid = move_left(grid)
        elif move == 'u':
            grid = move_up(grid)
        elif move == 'd':
            grid = move_down(grid)
        elif move == 'r':
            grid = move_right(grid)
    return grid

initial_grid = [
    [0, 32, 2, 2],
    [0, 0, 0, 0],
    [32, 0, 0, 0],
    [4, 4, 0, 0]
]

moves = 'ludl'
final_grid = apply_moves(initial_grid, moves)
print(final_grid)