def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def get_minor_diagonal(grid):
    return [grid[i][6-i] for i in range(7)]

def find_diagonal_letter(grid):
    # Find the letter that appears in the minor diagonal
    diag = get_minor_diagonal(grid)
    for letter in 'abcdefg':
        # Check if this letter can be used for all empty diagonal positions
        can_use = True
        for i, cell in enumerate(diag):
            if cell and cell != letter:
                can_use = False
                break
        if can_use:
            return letter
    return None

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal
    if row + col == 6:
        diag_letter = find_diagonal_letter(grid)
        if num != diag_letter:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal cells
    for i in range(7):
        if not grid[i][6-i]:
            return i, 6-i
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid):
    # First, determine the diagonal letter
    diag_letter = find_diagonal_letter(grid)
    if not diag_letter:
        return False
    
    # Fill all empty diagonal cells with this letter
    for i in range(7):
        if not grid[i][6-i]:
            grid[i][6-i] = diag_letter
    
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
]

if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")