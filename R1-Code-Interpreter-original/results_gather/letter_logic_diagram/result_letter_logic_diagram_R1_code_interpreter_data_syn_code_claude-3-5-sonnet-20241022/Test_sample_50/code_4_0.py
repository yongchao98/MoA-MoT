def analyze_diagonal(grid):
    # Check what letters are already on the minor diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i]:
            diagonal_letters.add(grid[i][6-i])
    return list(diagonal_letters)[0] if diagonal_letters else None

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
    if row + col == 6 and letter != diag_letter:
        return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if i + (6-i) == 6 and not grid[i][6-i]:
            return i, 6-i
    # Then fill rest
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid, diag_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # If on minor diagonal, only try diagonal letter
    if row + col == 6:
        letters = [diag_letter]
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['c','','','','b','',''],
    ['','f','e','b','','g',''],
    ['','','b','','','',''],
    ['','b','a','','c','d','f'],
    ['b','a','','','','',''],
    ['','g','','','f','e',''],
    ['','c','','f','','','']
]

# Find what letter must be on the diagonal
diag_letter = analyze_diagonal(initial)
if not diag_letter:
    # Try each possible letter for diagonal
    for letter in 'abcdefg':
        # Make a copy of the initial grid
        test_grid = [row[:] for row in initial]
        # Check if this letter works with pre-filled cells
        valid = True
        for i in range(7):
            if test_grid[i][6-i]:
                if test_grid[i][6-i] != letter:
                    valid = False
                    break
        if valid:
            diag_letter = letter
            break

if diag_letter:
    # Pre-fill diagonal with the determined letter
    for i in range(7):
        if not initial[i][6-i]:
            initial[i][6-i] = diag_letter
    
    # Solve
    if solve(initial, diag_letter):
        # Print solution
        result = []
        for row in initial:
            result.append(','.join(row))
        print('\n'.join(result))
    else:
        print("No solution exists")
else:
    print("No valid diagonal letter found")