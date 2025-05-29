def is_valid(grid, row, col, letter, initial_grid):
    # Check if we're trying to change a pre-filled cell
    if initial_grid[row][col] != '':
        return letter == initial_grid[row][col]
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        # Must match the first letter we find on the diagonal in the initial grid
        for i in range(7):
            j = 6 - i
            if initial_grid[i][j] != '':
                return letter == initial_grid[i][j]
    
    return True

def solve_grid(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # If cell is pre-filled in initial grid, use that value
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve_grid(grid, initial_grid, row, col + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', '', '', '', 'd', 'g', 'c'],
    ['', 'a', '', 'd', '', 'c', ''],
    ['a', '', 'd', 'g', '', '', ''],
    ['', '', 'g', 'c', '', '', 'a'],
    ['d', '', 'c', 'e', '', '', 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', '', 'd', 'g']
]

grid = [[''] * 7 for _ in range(7)]

if solve_grid(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")