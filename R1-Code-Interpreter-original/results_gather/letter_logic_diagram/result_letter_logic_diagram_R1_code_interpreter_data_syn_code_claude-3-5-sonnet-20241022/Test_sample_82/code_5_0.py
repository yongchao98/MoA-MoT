def is_valid_diagonal(grid, letter):
    # Check if letter conflicts with any pre-filled diagonal positions
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != letter:
            return False
    # Check if using this letter on diagonal would make rows/columns impossible
    for i in range(7):
        if letter in [grid[i][j] for j in range(7) if j != 6-i and grid[i][j] != '']:
            return False
        if letter in [grid[j][6-i] for j in range(7) if j != i and grid[j][6-i] != '']:
            return False
    return True

def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    return True

def find_next_empty(grid):
    # First fill non-diagonal positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return i, j
    return None

def solve_puzzle():
    # Initial grid
    initial_grid = [
        ['f','','e','','','',''],
        ['','e','','d','','','f'],
        ['','','','','g','f','a'],
        ['','d','b','','','a',''],
        ['d','b','','','a','e',''],
        ['b','g','','a','e','c','d'],
        ['','f','a','','','d','']
    ]
    
    # Try each possible diagonal letter
    for diag_letter in 'abcdefg':
        # Create a working copy
        grid = [row[:] for row in initial_grid]
        
        # Check if this diagonal letter is possible
        if not is_valid_diagonal(grid, diag_letter):
            continue
            
        # Fill all diagonal positions with this letter
        for i in range(7):
            if grid[i][6-i] == '':
                grid[i][6-i] = diag_letter
                
        def solve(grid):
            empty = find_next_empty(grid)
            if not empty:
                return True
                
            row, col = empty
            used_in_row = set(grid[row])
            used_in_col = set(grid[i][col] for i in range(7))
            available = set('abcdefg') - used_in_row - used_in_col
            
            for letter in available:
                if is_valid_move(grid, row, col, letter):
                    grid[row][col] = letter
                    if solve(grid):
                        return True
                    grid[row][col] = ''
            return False
        
        # Try to solve with this diagonal letter
        if solve(grid):
            # Verify solution
            for row in grid:
                print(','.join(row))
            return True
            
    return False

solve_puzzle()