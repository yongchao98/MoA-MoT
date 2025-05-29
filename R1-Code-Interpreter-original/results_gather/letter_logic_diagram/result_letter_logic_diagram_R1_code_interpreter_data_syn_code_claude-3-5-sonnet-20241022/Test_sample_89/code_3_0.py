def verify_grid(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    return True

def matches_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

initial_grid = [
    ['', '', 'a', 'g', '', 'b', 'c'],
    ['', 'a', '', '', '', 'c', 'd'],
    ['a', 'g', 'e', 'b', '', '', ''],
    ['g', '', 'b', '', 'd', 'f', 'a'],
    ['e', 'b', '', '', 'f', '', 'g'],
    ['', '', 'd', 'f', 'a', '', ''],
    ['', 'd', 'f', 'a', 'g', 'e', 'b']
]

# First, determine what letter must be on the minor diagonal
diagonal_letters = set()
for i in range(7):
    if initial_grid[i][6-i] != '':
        diagonal_letters.add(initial_grid[i][6-i])
diagonal_letter = list(diagonal_letters)[0] if diagonal_letters else 'f'

# Create a working grid and fill in the minor diagonal
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]
        elif i + j == 6:  # minor diagonal
            grid[i][j] = diagonal_letter

def solve(grid, row, col):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return verify_grid(grid) and matches_initial(grid, initial_grid)
    
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    
    for letter in 'abcdefg':
        if letter not in used_in_row and letter not in used_in_col:
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

if solve(grid, 0, 0):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")