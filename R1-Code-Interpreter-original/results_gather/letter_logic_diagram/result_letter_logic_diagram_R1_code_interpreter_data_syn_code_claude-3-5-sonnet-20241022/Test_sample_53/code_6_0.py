def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def can_fill_with_diagonal(grid, diag_letter):
    # Make a copy of the grid
    test_grid = [row[:] for row in grid]
    
    # Fill all minor diagonal positions with diag_letter
    for i in range(7):
        if test_grid[i][6-i] == '':
            if not is_valid(test_grid, i, 6-i, diag_letter):
                return False
            test_grid[i][6-i] = diag_letter
    
    # Try to solve the rest
    return solve_remaining(test_grid)

def solve_remaining(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # Skip if position is on minor diagonal (already filled)
    if row + col == 6:
        return solve_remaining(grid)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', 'f', ''],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

# Try each letter for the minor diagonal
solution = None
for diag_letter in 'abcdefg':
    test_grid = [row[:] for row in initial_grid]
    # First check if this letter can be used for minor diagonal
    valid = True
    for i in range(7):
        if test_grid[i][6-i] != '' and test_grid[i][6-i] != diag_letter:
            valid = False
            break
    
    if valid and can_fill_with_diagonal(test_grid, diag_letter):
        solution = test_grid
        break

if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")