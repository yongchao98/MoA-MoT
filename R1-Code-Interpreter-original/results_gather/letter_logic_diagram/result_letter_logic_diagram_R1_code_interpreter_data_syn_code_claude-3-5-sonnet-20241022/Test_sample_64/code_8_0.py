def verify_constraints(grid, initial_grid):
    # Verify initial values
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Verify rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

def is_valid_partial(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, initial_grid, pos=0):
    if pos == 49:
        return verify_constraints(grid, initial_grid)
    
    row = pos // 7
    col = pos % 7
    
    # If position is pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, pos + 1)
    
    # If on minor diagonal, use the letter from another diagonal position if exists
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter:
            if is_valid_partial(grid, row, col, diagonal_letter):
                grid[row][col] = diagonal_letter
                if solve(grid, initial_grid, pos + 1):
                    return True
            grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_partial(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")