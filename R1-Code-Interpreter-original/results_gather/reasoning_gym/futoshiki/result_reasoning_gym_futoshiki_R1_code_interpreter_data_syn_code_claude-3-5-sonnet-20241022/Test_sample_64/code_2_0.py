def verify_constraints(grid, row, col, num, h_constraints, v_constraints):
    # Row and column uniqueness
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check all horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '<':
            if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
                return False
            if grid[row][col-1] == 0:  # Check if any number could satisfy
                possible = False
                for n in range(1, 7):
                    if n < num and n not in [grid[row][i] for i in range(6)]:
                        possible = True
                        break
                if not possible:
                    return False
        elif h_constraints[row][col-1] == '>':
            if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
                return False
            if grid[row][col-1] == 0:
                possible = False
                for n in range(1, 7):
                    if n > num and n not in [grid[row][i] for i in range(6)]:
                        possible = True
                        break
                if not possible:
                    return False

    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row < 5 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 5 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False

    return True

def solve(grid, h_constraints, v_constraints):
    empty = None
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 7):
        if verify_constraints(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial setup with given numbers and constraints
grid = [
    [0, 0, 0, 2, 4, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 4, 0, 5, 0, 0],
    [0, 6, 3, 0, 0, 0],
    [0, 0, 6, 0, 3, 5],
    [0, 0, 5, 0, 0, 4]
]

h_constraints = [
    ['<', ' ', ' ', '<', ' '],
    [' ', ' ', ' ', '>', ' '],
    [' ', ' ', ' ', ' ', ' '],
    [' ', ' ', '>', '<', ' '],
    [' ', ' ', ' ', ' ', ' '],
    ['>', ' ', ' ', ' ', ' ']
]

v_constraints = [
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', '∨', ' ', '∨', ' ', ' ']
]

if solve(grid, h_constraints, v_constraints):
    for row in grid:
        print(' '.join(map(str, row)))