def verify_constraints(grid, initial_grid):
    # Check pre-filled positions
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def get_next_cell(grid, initial_grid):
    # Priority 1: Pre-filled positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial_grid[i][j] != '':
                return (i, j, [initial_grid[i][j]])
    
    # Priority 2: Minor diagonal
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    
    if diag_letter:
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i, [diag_letter])
    
    # Priority 3: First empty cell in minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            used_letters = set()
            for r in range(7):
                for c in range(7):
                    if grid[r][c] != '':
                        used_letters.add(grid[r][c])
            available = [l for l in 'abcdefg' if l not in used_letters]
            return (i, 6-i, available)
    
    # Priority 4: Regular cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, list('abcdefg'))
    
    return None

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
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, initial_grid):
    next_cell = get_next_cell(grid, initial_grid)
    if not next_cell:
        return verify_constraints(grid, initial_grid)
    
    row, col, possible_letters = next_cell
    
    for letter in possible_letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', 'g', 'b', '', '', '', ''],
    ['', '', '', '', 'c', '', ''],
    ['b', '', '', '', '', 'a', ''],
    ['f', 'd', 'c', '', '', '', ''],
    ['', 'c', 'e', 'a', '', 'b', ''],
    ['', '', '', 'g', '', 'f', 'd'],
    ['e', 'a', 'g', '', 'f', '', '']
]

grid = [[''] * 7 for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

if solve(grid, initial_grid):
    print("<<<")
    for row in grid:
        print(",".join(row))
    print(">>>")
else:
    print("No solution exists")