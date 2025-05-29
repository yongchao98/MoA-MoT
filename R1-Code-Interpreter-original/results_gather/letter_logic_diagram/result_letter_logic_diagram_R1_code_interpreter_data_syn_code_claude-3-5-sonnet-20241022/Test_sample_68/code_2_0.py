def find_diagonal_letter(grid):
    # Check existing letters on minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
    return None

def is_valid(grid, row, col, letter, diagonal_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # Position is on minor diagonal
        return letter == diagonal_letter
    
    return True

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    if row + col == 6:  # On minor diagonal
        letters = [diagonal_letter]
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# Find what letter must be on the diagonal
diagonal_letter = find_diagonal_letter(grid)
if diagonal_letter is None:
    print("No diagonal letter found in initial grid")
else:
    if solve_grid(grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution exists")