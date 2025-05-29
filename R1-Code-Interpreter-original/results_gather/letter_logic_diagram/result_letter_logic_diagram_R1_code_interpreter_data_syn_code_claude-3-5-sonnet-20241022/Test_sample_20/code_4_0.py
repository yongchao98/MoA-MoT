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

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # Skip if this position is on the minor diagonal (already filled with 'e')
    if row + col == 6:
        return solve(grid)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['b', '', '', 'g', '', 'c', 'e'],
    ['', 'd', 'g', 'a', '', 'e', ''],
    ['d', '', 'a', 'c', '', 'b', 'f'],
    ['g', '', '', 'e', '', '', 'd'],
    ['', '', 'e', 'b', '', '', ''],
    ['', 'e', '', 'f', 'd', 'g', 'a'],
    ['e', 'b', 'f', 'd', 'g', '', 'c']
]

# First, force all minor diagonal positions to be 'e'
for i in range(7):
    j = 6-i
    if grid[i][j] == '':
        grid[i][j] = 'e'
    elif grid[i][j] != 'e':
        print("Impossible: diagonal position already contains non-'e' value")
        exit()

if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")