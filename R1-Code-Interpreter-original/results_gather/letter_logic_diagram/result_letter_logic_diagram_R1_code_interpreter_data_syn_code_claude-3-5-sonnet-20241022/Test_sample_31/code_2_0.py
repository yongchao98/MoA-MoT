def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial_grid):
    # First check: if this position was pre-filled
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this cell is part of it
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve_grid(grid, initial_grid):
    # Find empty cell with most constraints first
    min_options = float('inf')
    empty = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Count available options for this cell
                options = sum(1 for letter in 'abcdefg' 
                            if is_valid(grid, i, j, letter, initial_grid))
                if options < min_options:
                    min_options = options
                    empty = (i, j)
    
    if not empty:
        return True
    
    row, col = empty
    letters = list('abcdefg')
    
    # If cell is on minor diagonal, prioritize existing diagonal letter
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = [grid[i][j]] + [l for l in letters if l != grid[i][j]]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

def verify_solution(grid, initial_grid):
    # Verify all initial constraints
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Verify rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7:  # Check row
            return False
        if len(set(grid[j][i] for j in range(7))) != 7:  # Check column
            return False
    
    # Verify minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

# Initial grid
initial_grid = [
    ['e', 'a', '', '', '', 'c', 'f'],
    ['', 'b', '', '', 'c', 'f', ''],
    ['b', 'g', 'd', 'c', 'f', '', 'a'],
    ['g', 'd', '', 'f', 'e', '', 'b'],
    ['d', '', '', 'e', '', '', 'g'],
    ['', 'f', 'e', '', 'b', '', 'd'],
    ['', '', 'a', '', 'g', 'd', 'c']
]

# Create working grid
grid = [row[:] for row in initial_grid]

if solve_grid(grid, initial_grid):
    if verify_solution(grid, initial_grid):
        print_grid(grid)
    else:
        print("Solution found but verification failed")
else:
    print("No solution exists")