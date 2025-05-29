def get_diagonal_candidates(grid):
    used_in_rows = set()
    used_in_cols = set()
    diagonal_positions = []
    
    # Check what letters are already used in rows and columns that intersect with diagonal
    for i in range(7):
        j = 6-i
        diagonal_positions.append((i,j))
        if grid[i][j] != '':
            return {grid[i][j]}  # If we have a diagonal letter, that's our only option
        
        # Add all letters in this row and column
        used_in_rows.update(x for x in grid[i] if x != '')
        used_in_cols.update(grid[r][j] for r in range(7) if grid[r][j] != '')
    
    # Potential diagonal letters are those not conflicting with any row or column
    all_letters = set('abcdefg')
    return all_letters - (used_in_rows | used_in_cols)

def create_grid_with_diagonal(original_grid, diagonal_letter):
    grid = [row[:] for row in original_grid]
    # Fill diagonal with the chosen letter
    for i in range(7):
        j = 6-i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
        elif grid[i][j] != diagonal_letter:
            return None  # Invalid diagonal letter
    return grid

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

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        # Skip if this would violate diagonal constraint
        if row + col == 6 and letter != grid[0][6]:  # Compare with first diagonal letter
            continue
            
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['a','e','f','g','','c','d'],
    ['e','f','g','b','','d','a'],
    ['','g','b','','','','e'],
    ['','b','','','a','','f'],
    ['','c','d','','e','','g'],
    ['c','d','a','e','f','g','b'],
    ['d','a','e','f','','b','']
]

# Get possible diagonal letters
diagonal_candidates = get_diagonal_candidates(initial_grid)

# Try each possible diagonal letter
for diagonal_letter in diagonal_candidates:
    grid = create_grid_with_diagonal(initial_grid, diagonal_letter)
    if grid is None:
        continue
        
    if solve_grid(grid):
        # Verify solution
        for i in range(7):
            if grid[i][6-i] != diagonal_letter:
                continue
            
        # Print solution
        for row in grid:
            print(','.join(row))
        break