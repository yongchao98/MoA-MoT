def find_diagonal_letter(grid):
    # Check pre-filled cells on minor diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    if len(diagonal_letters) > 1:
        return None  # Conflict in pre-filled diagonal
    elif len(diagonal_letters) == 1:
        return list(diagonal_letters)[0]
    return 'e'  # Default if no pre-filled diagonal cells

def is_valid(grid):
    # Check rows
    for row in grid:
        if '' in row:
            continue
        if len(set(row)) != 7:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if '' in col:
            continue
        if len(set(col)) != 7:
            return False
    
    # Check diagonal
    diag = [grid[i][6-i] for i in range(7)]
    diag_letters = set(x for x in diag if x != '')
    if len(diag_letters) > 1:
        return False
    
    return True

def solve(grid, diag_letter):
    if not is_valid(grid):
        return False

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
    
    i, j = empty
    # If on diagonal, must use diagonal letter
    if i + j == 6:
        if is_valid(grid):
            grid[i][j] = diag_letter
            if solve(grid, diag_letter):
                return True
        grid[i][j] = ''
        return False
    
    # Try each letter
    used_row = set(grid[i])
    used_col = set(grid[r][j] for r in range(7))
    
    for letter in 'abcdefg':
        if letter not in used_row and letter not in used_col:
            grid[i][j] = letter
            if solve(grid, diag_letter):
                return True
            grid[i][j] = ''
    
    return False

# Initialize grid
grid = [
    ['g','b','f','','c','a',''],
    ['','f','d','','a','e',''],
    ['','d','c','','','',''],
    ['','','','','','',''],
    ['','','e','','b','',''],
    ['a','e','','','','d','c'],
    ['e','','','','','c','']
]

# Find diagonal letter
diag_letter = find_diagonal_letter(grid)

# Fill diagonal first
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = diag_letter

if solve(grid, diag_letter):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")