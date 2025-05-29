def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_fixed_cells(grid, diag_letter):
    # Check if diagonal letter conflicts with pre-filled cells
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diag_letter:
            return False
        # Check if diagonal letter appears elsewhere in rows/columns
        for j in range(7):
            if j != 6-i:  # Skip diagonal position
                if grid[i][j] == diag_letter:  # Check row
                    return False
                if grid[j][i] == diag_letter:  # Check column
                    return False
    return True

def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
    return True

def solve(grid, diag_letter):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                if i + j == 6:  # Diagonal position
                    grid[i][j] = diag_letter
                    continue
                
                for letter in 'abcdefg':
                    if letter != diag_letter and is_valid(grid, i, j, letter):
                        grid[i][j] = letter
                        if solve(grid, diag_letter):
                            return True
                        grid[i][j] = ''
                return False
    return True

# Initial grid
initial = [
    ['', 'g', 'e', '', 'b', '', ''],
    ['g', 'e', 'f', 'b', '', '', ''],
    ['', '', '', '', '', '', 'g'],
    ['f', '', 'd', '', '', 'g', 'e'],
    ['', 'd', '', '', '', '', ''],
    ['d', 'c', '', '', '', '', 'b'],
    ['c', '', 'g', '', '', 'b', '']
]

# Try each letter as the diagonal letter
for diag_letter in 'abcdefg':
    # Create a copy of the initial grid
    test_grid = [row[:] for row in initial]
    
    # Check if this diagonal letter is possible with fixed cells
    if not check_fixed_cells(test_grid, diag_letter):
        continue
    
    # Fill all diagonal positions first
    for i in range(7):
        if test_grid[i][6-i] == '':
            test_grid[i][6-i] = diag_letter
    
    # Try to solve the rest of the grid
    if solve(test_grid, diag_letter):
        # Verify solution
        rows_valid = all(len(set(row)) == 7 for row in test_grid)
        cols_valid = all(len(set(col)) == 7 for col in zip(*test_grid))
        diag_valid = len(set(test_grid[i][6-i] for i in range(7))) == 1
        
        if rows_valid and cols_valid and diag_valid:
            print_grid(test_grid)
            break