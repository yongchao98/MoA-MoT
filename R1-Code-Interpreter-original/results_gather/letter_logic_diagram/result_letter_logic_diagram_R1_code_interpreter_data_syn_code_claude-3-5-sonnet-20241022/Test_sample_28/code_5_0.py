def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def is_valid(grid, row, col, num):
    # If it's a diagonal position, must be 'd'
    if row + col == 6 and num != 'd':
        return False
        
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    return True

def find_empty(grid):
    # First fill diagonal cells
    for i in range(7):
        if not grid[i][6-i]:
            return i, 6-i
    # Then other cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # If it's a diagonal position, only try 'd'
    letters = 'd' if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
]

# First fill all diagonal positions with 'd'
for i in range(7):
    if not grid[i][6-i]:
        grid[i][6-i] = 'd'

if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")