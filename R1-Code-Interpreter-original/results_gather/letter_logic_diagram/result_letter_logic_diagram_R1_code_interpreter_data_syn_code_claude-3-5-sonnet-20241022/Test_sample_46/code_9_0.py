def analyze_diagonal_constraints(grid):
    # Track what letters can't be on diagonal due to row/column conflicts
    cant_be_diagonal = set()
    
    # Check each diagonal position
    for i in range(7):
        # If position is filled, that must be our diagonal letter
        if grid[i][6-i] != '':
            return grid[i][6-i]
        
        # Check what letters appear in same row/column as this diagonal position
        for j in range(7):
            if grid[i][j] != '':  # Check row
                cant_be_diagonal.add(grid[i][j])
            if grid[j][6-i] != '':  # Check column
                cant_be_diagonal.add(grid[j][6-i])
    
    # Return possible diagonal letters
    possible = set('abcdefg') - cant_be_diagonal
    return list(possible)[0] if possible else None

def is_valid(grid, row, col, letter, diag_letter):
    # If on diagonal, must be diagonal letter
    if row + col == 6:
        return letter == diag_letter
    
    # Can't use diagonal letter off diagonal
    if letter == diag_letter:
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
    
    if not empty:
        return True
    
    row, col = empty
    
    # If diagonal position, only one choice
    if row + col == 6:
        if is_valid(grid, row, col, diag_letter, diag_letter):
            grid[row][col] = diag_letter
            return solve(grid, diag_letter)
        return False
    
    # Try each letter except diagonal letter
    letters = [l for l in 'abcdefg' if l != diag_letter]
    for letter in letters:
        if is_valid(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve(grid, diag_letter):
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

# First, determine what letter must be on diagonal
diag_letter = analyze_diagonal_constraints(initial)

if diag_letter:
    # Create working copy
    grid = [row[:] for row in initial]
    
    # Fill all diagonal positions first
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diag_letter
    
    # Solve rest of grid
    if solve(grid, diag_letter):
        # Print solution
        for row in grid:
            print(','.join(row))
    else:
        print("No solution exists")
else:
    print("Could not determine diagonal letter")