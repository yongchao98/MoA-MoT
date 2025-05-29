def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must match diagonal letter
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # If this is first empty minor diagonal position, try all letters
    # Otherwise, use the established minor diagonal letter
    if row + col == 6:
        minor_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                minor_letter = grid[i][6-i]
                break
        letters = 'abcdefg' if minor_letter is None else minor_letter
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
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', 'f', ''],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")