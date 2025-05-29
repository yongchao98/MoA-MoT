def get_initial_constraints(initial):
    constraints = []
    # Get all pre-filled positions
    for i in range(7):
        for j in range(7):
            if initial[i][j]:
                constraints.append((i, j, initial[i][j]))
    return constraints

def get_minor_diagonal():
    return [(i, 6-i) for i in range(7)]

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6-i
            if grid[i][j] and grid[i][j] != letter:
                return False
                
    return True

def solve(grid, pos, initial_constraints, diag_positions):
    if pos >= 49:
        return True
        
    row = pos // 7
    col = pos % 7
    
    # Skip if this position is pre-filled
    if any(r == row and c == col for r, c, _ in initial_constraints):
        return solve(grid, pos + 1, initial_constraints, diag_positions)
    
    # Get the set of letters to try
    if (row, col) in diag_positions:
        # If on minor diagonal, find existing diagonal letter
        diag_letter = None
        for i, j in diag_positions:
            if grid[i][j]:
                diag_letter = grid[i][j]
                break
        letters = [diag_letter] if diag_letter else list('abcdefg')
    else:
        letters = list('abcdefg')
    
    # Try each possible letter
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, pos + 1, initial_constraints, diag_positions):
                return True
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

# Get initial constraints and diagonal positions
initial_constraints = get_initial_constraints(initial)
diag_positions = get_minor_diagonal()

# Fill in initial constraints first
for row, col, letter in initial_constraints:
    grid[row][col] = letter

# Try to solve
if solve(grid, 0, initial_constraints, diag_positions):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")