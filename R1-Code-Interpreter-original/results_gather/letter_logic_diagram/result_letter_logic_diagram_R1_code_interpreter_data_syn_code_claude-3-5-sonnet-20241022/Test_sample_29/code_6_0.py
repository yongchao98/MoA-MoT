def analyze_diagonal_constraints(grid):
    # First check if any diagonal position is filled
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letter = grid[i][6-i]
            break
    
    if diagonal_letter:
        return diagonal_letter
    
    # If no diagonal position is filled, analyze what letters are possible
    possible = set('abcdefg')
    
    # Remove letters that can't be used due to row/column constraints
    for i in range(7):
        j = 6-i  # diagonal position
        # Check row constraints
        row_used = set(grid[i]) - {''}
        possible -= row_used
        # Check column constraints
        col_used = set(grid[r][j] for r in range(7)) - {''}
        possible -= col_used
    
    return list(possible)[0] if len(possible) == 1 else None

def create_base_grid(original_grid, diagonal_letter):
    grid = [row[:] for row in original_grid]
    # Fill all diagonal positions with diagonal_letter
    for i in range(7):
        j = 6-i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
        elif grid[i][j] != diagonal_letter:
            return None
    return grid

def is_valid_partial(grid):
    # Check rows
    for row in grid:
        seen = set()
        for cell in row:
            if cell != '' and cell in seen:
                return False
            seen.add(cell)
    
    # Check columns
    for j in range(7):
        seen = set()
        for i in range(7):
            if grid[i][j] != '' and grid[i][j] in seen:
                return False
            seen.add(grid[i][j])
    
    # Check diagonal
    diagonal = set(grid[i][6-i] for i in range(7) if grid[i][6-i] != '')
    if len(diagonal) > 1:
        return False
    
    return True

def solve_with_diagonal(grid):
    if not is_valid_partial(grid):
        return False
        
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i,j)
                break
        if empty:
            break
            
    if not empty:
        return True
        
    i, j = empty
    diagonal_letter = grid[0][6]  # Get diagonal letter from first diagonal position
    
    # Try each possible letter
    for letter in 'abcdefg':
        if i + j == 6 and letter != diagonal_letter:
            continue
            
        grid[i][j] = letter
        if is_valid_partial(grid) and solve_with_diagonal(grid):
            return True
        grid[i][j] = ''
    
    return False

# Initial grid
grid = [
    ['a','e','f','g','','c','d'],
    ['e','f','g','b','','d','a'],
    ['','g','b','','','','e'],
    ['','b','','','a','','f'],
    ['','c','d','','e','','g'],
    ['c','d','a','e','f','g','b'],
    ['d','a','e','f','','b','']
]

# Find the required diagonal letter
diagonal_letter = analyze_diagonal_constraints(grid)
if diagonal_letter:
    base_grid = create_base_grid(grid, diagonal_letter)
    if base_grid and solve_with_diagonal(base_grid):
        for row in base_grid:
            print(','.join(row))
    else:
        print("No solution exists")
else:
    print("Could not determine diagonal letter")