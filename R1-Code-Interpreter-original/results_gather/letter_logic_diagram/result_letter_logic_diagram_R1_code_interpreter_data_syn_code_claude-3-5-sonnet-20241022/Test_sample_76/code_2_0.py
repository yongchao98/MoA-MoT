def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:
        # Get the first letter on minor diagonal
        first_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                first_letter = grid[i][6-i]
                break
        if first_letter and letter != first_letter:
            return False
    
    return True

def print_grid(grid):
    for row in grid:
        print(','.join(row))

def solve_grid(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve_grid(grid, row, col + 1)

    # Try each letter
    letters = 'abcdefg'
    # If on minor diagonal, try only the first letter found on it
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break

    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with exact initial values
initial_grid = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

# Solve
print('<<<')
if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")
print('>>>')