def print_solution(grid):
    print("<<<")
    for row in grid:
        print(",".join(row))
    print(">>>")

def is_valid(grid, row, col, letter, initial_grid):
    # Must match pre-filled value
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
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
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid, initial_grid):
    # First priority: pre-filled positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial_grid[i][j] != '':
                return (i, j)
    
    # Second priority: diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Last: remaining positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    
    return None

def solve(grid, initial_grid):
    pos = find_empty(grid, initial_grid)
    if not pos:
        return True
    
    row, col = pos
    
    # If this is a pre-filled position, only try that letter
    if initial_grid[row][col] != '':
        letter = initial_grid[row][col]
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False
    
    # If this is on diagonal, try to find existing diagonal letter
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diag_letter = grid[i][j]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter, initial_grid):
                grid[row][col] = diag_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    letters = 'abcdefg'
    if row + col == 6 and not any(grid[i][6-i] != '' for i in range(7)):
        # If this is first diagonal position, prefer letters from pre-filled positions
        letters = [l for l in 'abcdefg' if any(l == initial_grid[i][j] 
                  for i in range(7) for j in range(7) if initial_grid[i][j] != '')]
    
    for letter in letters:
        if is_valid(grid, row, col, letter, initial_grid):
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
    print_solution(grid)
else:
    print("No solution exists")