def is_valid(grid, row, col, letter, initial_grid):
    # Must match initial grid if cell was pre-filled
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row - no duplicates
    temp_grid = [x for x in grid[row]]
    temp_grid[col] = letter
    if temp_grid.count(letter) > 1:
        return False
    
    # Check column - no duplicates
    col_values = [grid[i][col] for i in range(7)]
    col_values[row] = letter
    if col_values.count(letter) > 1:
        return False
    
    # Check minor diagonal
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
            if i == row and j == col:
                continue
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        # Verify final solution
        for i in range(7):
            row_letters = [x for x in grid[i] if x != '']
            if len(row_letters) != len(set(row_letters)):
                return False
            col_letters = [grid[j][i] for j in range(7) if grid[j][i] != '']
            if len(col_letters) != len(set(col_letters)):
                return False
        return True
    
    row, col = empty
    
    # If on minor diagonal, try to match existing diagonal letter
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter:
            if is_valid(grid, row, col, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'g', 'a', '', 'd', '', ''],
    ['', 'a', '', '', '', 'c', 'f'],
    ['a', '', 'd', 'b', 'c', '', 'g'],
    ['e', 'd', '', 'c', 'f', 'g', 'a'],
    ['', 'b', 'c', '', '', 'a', ''],
    ['', 'c', '', 'g', 'a', 'e', 'd'],
    ['', 'f', '', 'a', '', 'd', 'b']
]

# Create working copy
grid = [row[:] for row in initial_grid]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")