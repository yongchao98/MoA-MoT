def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # For minor diagonal cells, must be 'd'
    if row + col == 6 and letter != 'd':
        return False
    
    return True

def solve_grid(grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
    
    row = pos // 7
    col = pos % 7
    
    # If cell is pre-filled, move to next cell
    if grid[row][col] != '':
        return solve_grid(grid, pos + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with the given puzzle
initial_grid = [
    ['', 'e', '', '', '', 'g', ''],
    ['e', 'c', 'a', 'b', '', 'd', 'f'],
    ['', '', 'b', 'g', '', 'f', 'e'],
    ['', 'b', 'g', '', '', 'e', 'c'],
    ['b', '', 'd', 'f', '', 'c', ''],
    ['g', 'd', '', '', 'c', 'a', 'b'],
    ['d', '', 'e', '', 'a', 'b', '']
]

# First, set 'd' on empty minor diagonal cells
for i in range(7):
    if i + (6-i) == 6:  # on minor diagonal
        if initial_grid[i][6-i] == '':
            initial_grid[i][6-i] = 'd'

if solve_grid(initial_grid):
    print_grid(initial_grid)
else:
    print("No solution found")