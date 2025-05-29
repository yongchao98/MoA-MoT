def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid):
    # Check rows
    for row in grid:
        if row.count('') > 0:  # Skip incomplete rows
            continue
        if len(set(row)) != 7:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if col.count('') > 0:  # Skip incomplete columns
            continue
        if len(set(col)) != 7:
            return False
    
    # Check diagonal
    for i in range(7):
        if grid[i][6-i] != 'g':
            return False
    
    return True

def solve(grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return is_valid(grid)
    
    row = pos // 7
    col = pos % 7
    
    # If this position is pre-filled or on diagonal (which should be 'g')
    if grid[row][col] != '' or (row + col == 6):
        return solve(grid, pos + 1)
    
    for letter in 'abcdefg':
        # Try placing the letter
        grid[row][col] = letter
        if is_valid(grid):  # Check if current state is valid
            if solve(grid, pos + 1):
                return True
        grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['f', 'd', '', '', '', 'c', 'g'],
    ['d', 'e', '', 'a', 'c', '', 'f'],
    ['e', 'b', '', '', 'g', 'f', ''],
    ['', '', 'c', 'g', 'f', '', 'e'],
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['', 'g', '', 'd', 'e', '', 'a'],
    ['g', 'f', '', 'e', 'b', 'a', 'c']
]

# First, force 'g' on all diagonal positions
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'g'

# Now solve the rest
if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")