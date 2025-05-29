def print_grid(grid):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")

# Initialize the grid with pre-filled values
initial_grid = [
    ['', '', 'a', 'g', '', 'b', 'c'],
    ['', 'a', '', '', '', 'c', 'd'],
    ['a', 'g', 'e', 'b', '', '', ''],
    ['g', '', 'b', '', 'd', 'f', 'a'],
    ['e', 'b', '', '', 'f', '', 'g'],
    ['', '', 'd', 'f', 'a', '', ''],
    ['', 'd', 'f', 'a', 'g', 'e', 'b']
]

def is_valid_move(grid, row, col, num, initial_grid):
    # Check if this position was pre-filled
    if initial_grid[row][col] != '':
        return num == initial_grid[row][col]
    
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
            
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != num:
                return False
    
    return True

def find_empty_location(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_sudoku(grid, initial_grid):
    # Find empty location
    find = find_empty_location(grid)
    if not find:
        return True
    
    row, col = find
    
    # If this is on minor diagonal, we need to find what letter should be there
    if row + col == 6:
        for i in range(7):
            if initial_grid[i][6-i] != '':
                # Try only this letter
                if is_valid_move(grid, row, col, initial_grid[i][6-i], initial_grid):
                    grid[row][col] = initial_grid[i][6-i]
                    if solve_sudoku(grid, initial_grid):
                        return True
                    grid[row][col] = ''
                return False
    
    # Try digits 1 to 9
    for num in 'abcdefg':
        if is_valid_move(grid, row, col, num, initial_grid):
            grid[row][col] = num
            if solve_sudoku(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Create a working copy of the grid
grid = [row[:] for row in initial_grid]

# Solve the puzzle
if solve_sudoku(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")