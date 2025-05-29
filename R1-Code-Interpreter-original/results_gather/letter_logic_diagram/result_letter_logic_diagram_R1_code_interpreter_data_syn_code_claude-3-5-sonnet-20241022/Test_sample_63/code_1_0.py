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
    if row + col == 6:
        # Get the letter that should be on minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve_grid(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # If this position was pre-filled, skip it
    if initial_grid[row][col] != '':
        if grid[row][col] != initial_grid[row][col]:
            grid[row][col] = initial_grid[row][col]
        return solve_grid(grid, initial_grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with pre-filled values
initial_grid = [
    ['', 'g', 'b', '', '', '', ''],
    ['', '', '', '', 'c', '', ''],
    ['b', '', '', '', '', 'a', ''],
    ['f', 'd', 'c', '', '', '', ''],
    ['', 'c', 'e', 'a', '', 'b', ''],
    ['', '', '', 'g', '', 'f', 'd'],
    ['e', 'a', 'g', '', 'f', '', '']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Solve the puzzle
if solve_grid(grid, initial_grid):
    # Print solution in required format
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")
else:
    print("No solution exists")