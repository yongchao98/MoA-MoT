def get_diagonal_letter(grid):
    # Get the letter that must be on the minor diagonal
    # by checking pre-filled cells on the diagonal
    for i in range(7):
        j = 6-i
        if grid[i][j]:
            return grid[i][j]
    return None

def is_valid(grid, row, col, letter, diag_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        if letter != diag_letter:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal cells
    for i in range(7):
        j = 6-i
        if not grid[i][j]:
            return i, j
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid):
    # First, determine what letter must be on the diagonal
    diag_letter = get_diagonal_letter(grid)
    if not diag_letter:
        # If no diagonal letter is pre-filled, try each possibility
        for letter in 'abcdefg':
            if solve_with_diagonal(grid, letter):
                return True
        return False
    return solve_with_diagonal(grid, diag_letter)

def solve_with_diagonal(grid, diag_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    if row + col == 6:  # If on minor diagonal
        if is_valid(grid, row, col, diag_letter, diag_letter):
            grid[row][col] = diag_letter
            if solve_with_diagonal(grid, diag_letter):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diag_letter):
                grid[row][col] = letter
                if solve_with_diagonal(grid, diag_letter):
                    return True
                grid[row][col] = ''
    return False

# Initial grid
grid = [
    ['a', 'g', '', '', 'f', '', 'b'],
    ['', 'd', '', '', 'e', '', 'a'],
    ['d', '', 'f', 'e', '', 'a', 'g'],
    ['c', 'f', 'e', 'b', 'a', 'g', ''],
    ['', '', 'b', '', 'g', 'd', 'c'],
    ['e', 'b', 'a', 'g', 'd', 'c', ''],
    ['b', '', 'g', 'd', 'c', 'f', 'e']
]

if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")