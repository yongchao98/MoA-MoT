def can_use_diagonal_letter(grid, letter):
    # Check if letter can be used in all minor diagonal positions
    for i in range(7):
        j = 6 - i
        # If cell is filled and different from our letter, it won't work
        if grid[i][j] != '' and grid[i][j] != letter:
            return False
        # Check row conflicts
        for k in range(7):
            if k != j and grid[i][k] == letter:
                return False
        # Check column conflicts
        for k in range(7):
            if k != i and grid[k][j] == letter:
                return False
    return True

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

def solve_grid(grid):
    # First find which letter can be used for the minor diagonal
    diagonal_letter = None
    for letter in 'abcdefg':
        if can_use_diagonal_letter(grid, letter):
            diagonal_letter = letter
            break
    
    if diagonal_letter is None:
        return False
    
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
    
    # Now solve the rest of the grid
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', '', 'e', 'd', 'c', ''],
    ['f', '', 'e', '', '', '', ''],
    ['', '', 'd', '', 'g', '', ''],
    ['e', '', 'c', '', 'b', 'f', ''],
    ['d', '', '', 'b', '', 'a', 'e'],
    ['c', 'g', 'b', '', 'a', 'e', 'd'],
    ['', 'b', 'f', '', '', 'd', '']
]

if solve_grid(grid):
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print("<<<")
    print(result.strip())
    print(">>>")
else:
    print("No solution exists")