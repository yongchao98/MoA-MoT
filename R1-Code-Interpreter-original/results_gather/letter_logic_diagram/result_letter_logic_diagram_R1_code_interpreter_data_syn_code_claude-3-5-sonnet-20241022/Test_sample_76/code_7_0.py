def print_solution(grid):
    for row in grid:
        print(','.join(row))

def get_minor_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Check minor diagonal
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False

    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    diag_letter = None
    
    # If on minor diagonal, get the required letter
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter:
            if is_valid(grid, row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve(grid):
                    return True
                grid[row][col] = ''
            return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with initial values
grid = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

# First verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != '':
            letter = grid[i][j]
            grid[i][j] = ''
            if not is_valid(grid, i, j, letter):
                valid = False
            grid[i][j] = letter
            if not valid:
                break
    if not valid:
        break

print('<<<')
if valid and solve(grid):
    print_solution(grid)
else:
    print("No valid solution exists")
print('>>>')