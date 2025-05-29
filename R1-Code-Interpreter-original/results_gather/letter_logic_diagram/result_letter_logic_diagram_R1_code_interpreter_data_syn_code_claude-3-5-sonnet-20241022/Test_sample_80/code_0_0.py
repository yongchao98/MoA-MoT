def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7) if grid[i][col] != '']:
        return False
    
    # Check minor diagonal if this cell is on it
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', '', '', '', '', 'a', ''],
    ['', 'e', '', 'd', '', '', ''],
    ['', 'c', '', 'a', '', 'g', 'f'],
    ['', '', '', 'b', '', 'f', 'e'],
    ['d', 'a', '', '', '', 'e', ''],
    ['', '', '', 'f', 'e', '', ''],
    ['b', 'g', 'f', 'e', '', '', '']
]

if solve(initial_grid):
    print_grid(initial_grid)
else:
    print("No solution found")