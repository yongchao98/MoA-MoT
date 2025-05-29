def is_valid(grid, row, col, letter):
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
        # Get the first non-empty cell on minor diagonal
        first_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                first_letter = grid[i][6-i]
                break
        if first_letter and letter != first_letter:
            return False
    
    return True

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
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
grid = [
    ['', '', '', '', 'd', 'g', 'c'],
    ['', 'a', '', 'd', '', 'c', ''],
    ['a', '', 'd', 'g', '', '', ''],
    ['', '', 'g', 'c', '', '', 'a'],
    ['d', '', 'c', 'e', '', '', 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', '', 'd', 'g']
]

if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")