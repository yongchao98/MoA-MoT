def is_valid_non_diagonal(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def fill_diagonal_first(grid):
    # First, fill all minor diagonal positions with 'c'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'c'
        elif grid[i][6-i] != 'c':
            return False
    return True

def solve_remaining(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip if cell is already filled or on minor diagonal
    if grid[row][col] != '' or row + col == 6:
        return solve_remaining(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if letter != 'c' and is_valid_non_diagonal(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', 'a', '', 'g', 'b', '', ''],
    ['a', 'd', 'g', '', '', 'c', ''],
    ['d', '', '', 'e', 'c', '', 'a'],
    ['', 'b', 'e', '', 'f', 'a', 'd'],
    ['b', 'e', 'c', '', '', 'd', 'g'],
    ['', '', '', '', 'd', 'g', ''],
    ['', '', 'a', '', '', '', '']
]

# Make a copy of the original grid
original_grid = [row[:] for row in grid]

# First phase: Fill diagonal with 'c'
if fill_diagonal_first(grid):
    # Second phase: Fill remaining cells
    if solve_remaining(grid):
        result = []
        for row in grid:
            result.append(','.join(row))
        print('<<<')
        for row in result:
        
            print(row)
        print('>>>')
    else:
        print("No solution exists after filling diagonal")
else:
    print("Cannot fill diagonal with 'c' due to conflicts")