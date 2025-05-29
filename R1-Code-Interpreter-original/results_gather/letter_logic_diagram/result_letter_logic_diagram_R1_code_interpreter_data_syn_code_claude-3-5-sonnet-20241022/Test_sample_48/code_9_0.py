def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    # First fill remaining diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i, True)  # True indicates diagonal position
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)  # False indicates non-diagonal position
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col, is_diagonal = pos
    
    if is_diagonal:
        # Only try 'g' for diagonal positions
        if is_valid_move(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve(grid):
                return True
            grid[row][col] = ''
    else:
        # Try all letters except those already in row/column
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['f', 'd', '', '', '', 'c', 'g'],
    ['d', 'e', '', 'a', 'c', '', 'f'],
    ['e', 'b', '', '', 'g', 'f', 'g'],  # Force g on diagonal
    ['', '', 'c', 'g', 'f', 'g', 'e'],  # Force g on diagonal
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['', 'g', '', 'd', 'e', 'g', 'a'],  # Force g on diagonal
    ['g', 'f', '', 'e', 'b', 'a', 'c']
]

# First ensure all diagonal positions have 'g'
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'g'

if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")