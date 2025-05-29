def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint if cell is on minor diagonal
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

def solve(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, i, j, letter):
                        grid[i][j] = letter
                        if solve(grid):
                            return True
                        grid[i][j] = ''
                return False
    return True

# Initialize grid
initial_grid = [
    ['f','','','c','a','','g'],
    ['e','b','c','a','d','g',''],
    ['','c','','','g','f','e'],
    ['c','a','','','','e','b'],
    ['a','d','g','','e','','c'],
    ['d','g','','e','','c',''],
    ['','','','b','','a','']
]

# Solve
if solve(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")