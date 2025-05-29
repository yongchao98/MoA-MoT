def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_diagonal_candidates(grid):
    # Get all possible letters that could work on the diagonal
    # based on pre-filled cells and available positions
    used_in_diag = set()
    for i in range(7):
        if grid[i][6-i]:
            used_in_diag.add(grid[i][6-i])
    
    if len(used_in_diag) > 1:
        return []  # Conflict in diagonal
    elif len(used_in_diag) == 1:
        return list(used_in_diag)
    
    return list('abcdefg')

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
            
    return True

def solve_with_diagonal(grid, diag_letter):
    def find_empty():
        for i in range(7):
            for j in range(7):
                if not grid[i][j]:
                    return i, j
        return None
    
    pos = find_empty()
    if not pos:
        return True
        
    row, col = pos
    
    # If this is a diagonal position, we must use diag_letter
    if row + col == 6:
        if is_valid_placement(grid, row, col, diag_letter):
            grid[row][col] = diag_letter
            if solve_with_diagonal(grid, diag_letter):
                return True
            grid[row][col] = ''
    else:
        # Try each letter except diag_letter if this position needs it for row/column completion
        for letter in 'abcdefg':
            if is_valid_placement(grid, row, col, letter):
                grid[row][col] = letter
                if solve_with_diagonal(grid, diag_letter):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['c','','','','b','',''],
    ['','f','e','b','','g',''],
    ['','','b','','','',''],
    ['','b','a','','c','d','f'],
    ['b','a','','','','',''],
    ['','g','','','f','e',''],
    ['','c','','f','','','']
]

# Try each possible diagonal letter
for diag_letter in get_diagonal_candidates(initial):
    # Make a copy of the initial grid
    test_grid = [row[:] for row in initial]
    
    # Pre-fill diagonal with the current letter
    valid = True
    for i in range(7):
        if test_grid[i][6-i]:
            if test_grid[i][6-i] != diag_letter:
                valid = False
                break
        else:
            test_grid[i][6-i] = diag_letter
    
    if not valid:
        continue
        
    if solve_with_diagonal(test_grid, diag_letter):
        # Print solution
        result = []
        for row in test_grid:
            result.append(','.join(row))
        print('\n'.join(result))
        break
else:
    print("No solution exists")