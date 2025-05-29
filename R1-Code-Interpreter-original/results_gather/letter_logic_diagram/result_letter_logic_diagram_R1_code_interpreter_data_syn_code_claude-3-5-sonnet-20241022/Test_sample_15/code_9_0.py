def get_possible_letters(grid, row, col, diagonal_letter):
    # Get all possible letters for this position
    used = set()
    
    # Check row
    for j in range(7):
        if grid[row][j] != '':
            used.add(grid[row][j])
    
    # Check column
    for i in range(7):
        if grid[i][col] != '':
            used.add(grid[i][col])
    
    # If on diagonal, must be diagonal_letter
    if row + col == 6:
        return {diagonal_letter} - used
    
    # Otherwise, can't be diagonal_letter
    possible = set('abcdefg') - used
    if row + col != 6:  # If not on diagonal
        possible.discard(diagonal_letter)  # Can't use diagonal letter
    
    return possible

def find_most_constrained_cell(grid):
    min_possibilities = float('inf')
    best_cell = None
    best_possible = None
    
    # First check non-diagonal positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                possible = get_possible_letters(grid, i, j, 'e')  # Try 'e' as diagonal letter
                if len(possible) < min_possibilities:
                    min_possibilities = len(possible)
                    best_cell = (i, j)
                    best_possible = possible
    
    return best_cell, best_possible

def is_valid_grid(grid):
    # Check rows
    for row in grid:
        if '' in row:
            continue
        if len(set(row)) != 7:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if '' in col:
            continue
        if len(set(col)) != 7:
            return False
    
    # Check diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if '' not in diagonal:
        if len(set(diagonal)) != 1:
            return False
    
    return True

def solve(grid):
    if not is_valid_grid(grid):
        return False
    
    cell, possible = find_most_constrained_cell(grid)
    if not cell:  # Grid is complete
        return True
    
    row, col = cell
    for letter in possible:
        grid[row][col] = letter
        if solve(grid):
            return True
        grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['c','','','','a','','b'],
    ['g','e','f','a','','','e'],
    ['','','','d','','e','g'],
    ['f','a','','e','c','g',''],
    ['','d','e','c','','e','f'],
    ['d','e','','','','f',''],
    ['e','','','','','','']
]

# Fill diagonal with 'e'
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'e'

if solve(grid):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")