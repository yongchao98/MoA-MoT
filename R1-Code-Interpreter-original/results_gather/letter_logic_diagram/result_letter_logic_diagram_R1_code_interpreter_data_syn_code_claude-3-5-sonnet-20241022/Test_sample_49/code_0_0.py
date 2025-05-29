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

def solve_grid(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    if grid[row][col] != '':
        return solve_grid(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', '', '', 'd', '', 'e'],
    ['c', 'f', '', '', '', 'g', ''],
    ['f', 'b', '', '', '', 'a', 'c'],
    ['', 'd', '', 'e', '', 'c', ''],
    ['', 'g', '', 'a', '', '', 'b'],
    ['', '', 'a', 'c', '', '', ''],
    ['', '', '', '', 'b', 'd', '']
]

# Solve the puzzle
if solve_grid(initial_grid):
    # Print solution in required format
    print('<<<')
    for row in initial_grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")