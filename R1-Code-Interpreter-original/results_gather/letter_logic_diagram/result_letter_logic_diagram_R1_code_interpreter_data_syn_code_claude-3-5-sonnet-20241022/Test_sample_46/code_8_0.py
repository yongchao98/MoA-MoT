def is_valid_diagonal(grid, letter):
    # Check if letter conflicts with pre-filled diagonal cells
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != letter:
            return False
        
    # Check if letter appears in non-diagonal positions
    for i in range(7):
        for j in range(7):
            if i + j != 6 and grid[i][j] == letter:
                return False
    return True

def is_valid_placement(grid, row, col, letter, diagonal_letter):
    # If this is a diagonal position, must be diagonal_letter
    if row + col == 6:
        return letter == diagonal_letter
    
    # Cannot use diagonal_letter in non-diagonal positions
    if letter == diagonal_letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    return True

def solve_grid(grid, diagonal_letter):
    empty = None
    # Find an empty cell (non-diagonal first)
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty non-diagonal cells, find empty diagonal cell
    if not empty:
        for i in range(7):
            if grid[i][6-i] == '':
                empty = (i, 6-i)
                break
    
    # If no empty cells, we're done
    if not empty:
        return True
        
    row, col = empty
    
    # If this is a diagonal position, only try diagonal_letter
    if row + col == 6:
        if is_valid_placement(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    else:
        # Try each letter except diagonal_letter
        for letter in 'abcdefg':
            if letter != diagonal_letter and is_valid_placement(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if solve_grid(grid, diagonal_letter):
                    return True
                grid[row][col] = ''
    
    return False

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

# Try each possible diagonal letter
for diagonal_letter in 'abcdefg':
    if not is_valid_diagonal(initial, diagonal_letter):
        continue
        
    # Create a working copy
    grid = [row[:] for row in initial]
    
    # Fill all empty diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    # Try to solve the rest
    if solve_grid(grid, diagonal_letter):
        # Verify solution
        rows_ok = all(len(set(row)) == 7 for row in grid)
        cols_ok = all(len(set(col)) == 7 for col in zip(*grid))
        diag_ok = all(grid[i][6-i] == diagonal_letter for i in range(7))
        
        if rows_ok and cols_ok and diag_ok:
            for row in grid:
                print(','.join(row))
            break