def validate_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] and grid[i][j] != initial[i][j]:
                return False
    return True

def get_next_cell(grid, initial):
    # First priority: Fill initial cells
    for i in range(7):
        for j in range(7):
            if initial[i][j] and not grid[i][j]:
                return (i, j)
    
    # Second priority: Fill minor diagonal
    for i in range(7):
        if not grid[i][6-i]:
            return (i, 6-i)
    
    # Third priority: Fill remaining cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return (i, j)
    
    return None

def get_candidates(grid, row, col, initial):
    if initial[row][col]:
        return {initial[row][col]}
    
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    candidates = set('abcdefg') - used_row - used_col
    
    # If on minor diagonal, must match existing diagonal letter
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                return {grid[i][6-i]} & candidates
    
    return candidates

def solve(grid, initial):
    cell = get_next_cell(grid, initial)
    if not cell:
        return validate_initial(grid, initial)
    
    row, col = cell
    candidates = get_candidates(grid, row, col, initial)
    
    # Try each candidate in a deterministic order
    for letter in sorted(candidates):
        # Place the letter
        grid[row][col] = letter
        
        # If this is on minor diagonal, propagate the letter
        if row + col == 6:
            # Save current state of diagonal
            old_diagonal = [grid[i][6-i] for i in range(7)]
            # Fill entire diagonal with this letter
            valid_diagonal = True
            for i in range(7):
                if not grid[i][6-i]:
                    grid[i][6-i] = letter
                elif grid[i][6-i] != letter:
                    valid_diagonal = False
                    break
            
            if valid_diagonal and solve(grid, initial):
                return True
                
            # Restore diagonal state
            for i in range(7):
                grid[i][6-i] = old_diagonal[i]
        else:
            if solve(grid, initial):
                return True
        
        # Backtrack
        grid[row][col] = ''
    
    return False

# Initialize puzzle
initial = [
    ['c','','g','','f','e',''],
    ['a','g','','','','d',''],
    ['','','','','d','',''],
    ['','','','d','','a','g'],
    ['','e','','c','a','g','b'],
    ['e','','','a','g','b','f'],
    ['d','','a','','b','','']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

if solve(grid, initial):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")