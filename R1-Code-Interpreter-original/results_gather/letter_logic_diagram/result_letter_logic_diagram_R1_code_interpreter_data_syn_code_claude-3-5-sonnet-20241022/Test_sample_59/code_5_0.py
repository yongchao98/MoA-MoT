def get_prefilled_positions(grid):
    prefilled = []
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                prefilled.append((i, j, grid[i][j]))
    return prefilled

def is_valid_placement(grid, row, col, letter, prefilled):
    # Check if this position was prefilled
    for r, c, l in prefilled:
        if r == row and c == col and l != letter:
            return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must be 'a'
    if row + col == 6 and letter != 'a':
        return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, prefilled):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    possible_letters = ['a'] if row + col == 6 else 'abcdefg'
    
    for letter in possible_letters:
        if is_valid_placement(grid, row, col, letter, prefilled):
            grid[row][col] = letter
            if solve(grid, prefilled):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'e', '', 'g', 'c', 'f', 'a'],
    ['', '', '', '', '', '', 'b'],
    ['', '', '', 'f', '', '', 'e'],
    ['', '', '', 'a', '', 'e', ''],
    ['', '', '', 'b', '', 'd', ''],
    ['f', '', 'b', 'e', '', '', 'c'],
    ['', '', '', 'd', 'g', '', 'f']
]

prefilled = get_prefilled_positions(initial_grid)

# First, verify that prefilled positions on minor diagonal are 'a'
for r, c, l in prefilled:
    if r + c == 6 and l != 'a':
        print("No solution exists - prefilled diagonal positions must be 'a'")
        exit()

if solve(initial_grid, prefilled):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")