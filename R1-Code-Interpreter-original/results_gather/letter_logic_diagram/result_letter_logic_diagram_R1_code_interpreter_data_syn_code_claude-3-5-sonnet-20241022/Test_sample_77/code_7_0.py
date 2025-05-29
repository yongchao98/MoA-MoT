def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_move(grid, row, col, letter, diag_letter):
    # If this is a minor diagonal position, must use diag_letter
    if row + col == 6:
        return letter == diag_letter
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve(grid, diag_letter):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty cell, we've found a solution
    if not empty:
        return True
    
    row, col = empty
    # If this is a minor diagonal position, only try diag_letter
    if row + col == 6:
        if is_valid_move(grid, row, col, diag_letter, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
    else:
        # Try all letters for non-diagonal positions
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter, diag_letter):
                grid[row][col] = letter
                if solve(grid, diag_letter):
                    return True
                grid[row][col] = ''
    
    return False

# Try each possible letter for the minor diagonal
initial_grid = [
    ['', '', 'd', 'c', 'a', '', ''],
    ['', 'd', 'c', 'a', '', 'f', ''],
    ['d', '', '', '', '', '', 'g'],
    ['', '', '', '', '', '', 'd'],
    ['', '', '', '', '', '', ''],
    ['b', 'f', '', 'g', 'd', '', 'a'],
    ['f', '', 'g', '', '', 'a', 'b']
]

for diag_letter in 'abcdefg':
    # Check if this letter conflicts with any pre-filled diagonal positions
    valid = True
    for i in range(7):
        if initial_grid[i][6-i] != '' and initial_grid[i][6-i] != diag_letter:
            valid = False
            break
    
    if not valid:
        continue
    
    # Create a fresh grid for this attempt
    grid = [row[:] for row in initial_grid]
    
    # Fill all empty minor diagonal positions with diag_letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diag_letter
    
    # Try to solve the rest of the grid
    if solve(grid, diag_letter):
        print_solution(grid)
        break