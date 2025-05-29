def verify_initial_constraints(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid(grid, row, col, val, initial):
    # Must match initial constraint if exists
    if initial[row][col] != '' and initial[row][col] != val:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == val:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == val:
            return False
    
    # Check minor diagonal - must be 'g' (from initial constraints)
    if row + col == 6 and val != 'g':
        return False
    
    return True

def solve(grid, initial):
    def find_empty():
        # First try minor diagonal
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i)
        # Then other positions
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    pos = find_empty()
    if not pos:
        return True
    
    row, col = pos
    candidates = ['g'] if row + col == 6 else 'abcdefg'
    
    for val in candidates:
        if is_valid(grid, row, col, val, initial):
            grid[row][col] = val
            if solve(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['b','','e','','','',''],
    ['','','','','f','',''],
    ['','','','f','g','',''],
    ['c','','f','g','b','d',''],
    ['','f','','b','d','','c'],
    ['','g','','d','','','a'],
    ['g','','','','','','f']
]

# Create working grid
grid = [row[:] for row in initial]

if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")