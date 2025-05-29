def get_minor_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If this is a minor diagonal position
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def find_empty(grid):
    # First try to fill minor diagonal
    minor_letter = get_minor_diagonal_letter(grid)
    if minor_letter:
        for i in range(7):
            if grid[i][6-i] == '' and is_valid(grid, i, 6-i, minor_letter):
                return i, 6-i
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, only try the existing diagonal letter
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        letters = [diag_letter] if diag_letter else list('abcdefg')
    else:
        letters = list('abcdefg')
    
    # Try 'f' first for minor diagonal (based on given constraints)
    if row + col == 6 and not diag_letter:
        letters = ['f'] + [l for l in letters if l != 'f']
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['d','','','e','b','','g'],
    ['','','e','b','','','d'],
    ['','','b','','','','c'],
    ['','','a','','d','c','f'],
    ['','a','g','d','','','e'],
    ['','','','','','',''],
    ['','','','f','','','a']
]

if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")