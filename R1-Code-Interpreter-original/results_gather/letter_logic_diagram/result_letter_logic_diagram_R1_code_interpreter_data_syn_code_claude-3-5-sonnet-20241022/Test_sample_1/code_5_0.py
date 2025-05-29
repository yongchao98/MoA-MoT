def is_valid(grid, row, col, letter, initial):
    # Must match initial if filled
    if initial[row][col] and initial[row][col] != letter:
        return False
    
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal - all must be the same
    if row + col == 6:
        diag = [(i, 6-i) for i in range(7)]
        for i, j in diag:
            if grid[i][j] and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill remaining positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # If on diagonal, try each letter and force it for all diagonal positions
    if row + col == 6:
        for letter in 'abcdefg':
            valid = True
            saved_values = []
            
            # Try to fill all diagonal positions with this letter
            for i in range(7):
                j = 6-i
                if grid[i][j] == '':
                    if is_valid(grid, i, j, letter, initial):
                        saved_values.append((i, j))
                        grid[i][j] = letter
                    else:
                        valid = False
                        break
                elif grid[i][j] != letter:
                    valid = False
                    break
            
            if valid and solve(grid, initial):
                return True
                
            # Backtrack
            for i, j in saved_values:
                grid[i][j] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial):
                grid[row][col] = letter
                if solve(grid, initial):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['b','','','c','g','f','a'],
    ['','e','c','g','f','a','b'],
    ['e','c','g','f','a','b','d'],
    ['c','','f','','b','','e'],
    ['g','f','','b','d','',''],
    ['','','b','','e','c','g'],
    ['','','d','e','','g','f']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

if solve(grid, initial):
    for row in grid:
        print(','.join(row))