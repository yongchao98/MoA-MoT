def get_diagonal_candidates(grid):
    # Get all letters that could potentially be on the diagonal
    # based on pre-filled positions
    candidates = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            candidates = {grid[i][6-i]}
            break
    return candidates

def is_valid_complete(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    return True

def is_valid_partial(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal consistency
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, diag_letter):
    def find_empty():
        # First fill diagonal positions
        for i in range(7):
            if grid[i][6-i] == '':
                return i, 6-i
        # Then other positions
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    pos = find_empty()
    if not pos:
        return is_valid_complete(grid)
    
    row, col = pos
    
    # If this is a diagonal position, only try the diagonal letter
    if row + col == 6:
        if is_valid_partial(grid, row, col, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid_partial(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, diag_letter):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['d','','','e','b','','g'],
    ['','','e','b','','','d'],
    ['','','b','','','','c'],
    ['','','a','','d','c','f'],
    ['','a','g','d','','','e'],
    ['','','','','','',''],
    ['','','','f','','','a']
]

# Try each possible diagonal letter
for diag_letter in get_diagonal_candidates(initial_grid):
    grid = [row[:] for row in initial_grid]
    
    # Pre-fill diagonal with the chosen letter where empty
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diag_letter
    
    if solve(grid, diag_letter):
        result = []
        for row in grid:
            result.append(','.join(row))
        print('\n'.join(result))
        break