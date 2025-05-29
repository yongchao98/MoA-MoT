def get_diagonal_letter(grid):
    # Find the first filled letter in the minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # If this is a minor diagonal position, must match diagonal letter
    if row + col == 6:
        diag_letter = get_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
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

def find_empty(grid):
    # First try to fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, only try the diagonal letter
    if row + col == 6:
        diag_letter = get_diagonal_letter(grid)
        letters = [diag_letter] if diag_letter else 'abcdefg'
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['c', 'g', '', 'd', 'e', 'b', ''],
    ['g', 'a', 'd', 'e', 'b', '', 'c'],
    ['', 'd', 'e', '', 'f', '', 'g'],
    ['d', 'e', '', 'f', '', '', ''],
    ['e', '', '', 'c', 'g', 'a', 'd'],
    ['b', 'f', 'c', 'g', 'a', '', ''],
    ['f', 'c', 'g', 'a', '', 'e', 'b']
]

# Find existing letter on minor diagonal
diag_letter = None
for i in range(7):
    if grid[i][6-i] != '':
        diag_letter = grid[i][6-i]
        break

if diag_letter is None:
    print("No diagonal letter found in initial grid")
else:
    if solve_grid(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution exists")