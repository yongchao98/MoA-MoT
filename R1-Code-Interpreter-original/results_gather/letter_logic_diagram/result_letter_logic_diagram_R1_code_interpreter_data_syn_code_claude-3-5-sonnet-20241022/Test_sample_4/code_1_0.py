def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # Position is on minor diagonal
        # Find the required diagonal letter
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If this is a diagonal position and we have a diagonal letter already
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['b', '', '', '', 'g', '', 'e'],
    ['', '', 'a', 'g', '', '', ''],
    ['f', 'a', '', '', '', 'b', 'd'],
    ['', 'g', '', '', '', '', ''],
    ['g', '', '', 'b', '', '', ''],
    ['', 'e', 'b', 'd', 'f', 'a', ''],
    ['', 'b', '', '', 'a', '', '']
]

if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")