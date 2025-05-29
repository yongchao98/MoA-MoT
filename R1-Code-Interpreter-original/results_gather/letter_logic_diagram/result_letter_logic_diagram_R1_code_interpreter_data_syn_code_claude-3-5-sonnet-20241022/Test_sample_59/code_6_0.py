def is_valid(grid, row, col, letter):
    # If this is a minor diagonal position, it must be 'a'
    if row + col == 6:
        return letter == 'a'
    
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
    # First fill minor diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid):
    # First, verify all pre-filled minor diagonal positions are 'a'
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != 'a':
            return False
    
    # Fill all empty minor diagonal positions with 'a'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'a'
    
    return solve_remaining(grid)

def solve_remaining(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid):
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

# Create a copy to work with
grid = [row[:] for row in initial_grid]

if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")