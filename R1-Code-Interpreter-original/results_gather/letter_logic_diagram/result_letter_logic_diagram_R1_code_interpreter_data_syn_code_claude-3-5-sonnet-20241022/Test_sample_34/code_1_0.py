def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve(grid, original):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # If it's on minor diagonal, use the first letter found on minor diagonal
                if i + j == 6:
                    diag_letter = None
                    for k in range(7):
                        if original[k][6-k] != '':
                            diag_letter = original[k][6-k]
                            break
                    if diag_letter:
                        if is_valid(grid, i, j, diag_letter):
                            grid[i][j] = diag_letter
                            if solve(grid, original):
                                return True
                            grid[i][j] = ''
                        return False
                
                for letter in 'abcdefg':
                    if is_valid(grid, i, j, letter):
                        grid[i][j] = letter
                        if solve(grid, original):
                            return True
                        grid[i][j] = ''
                return False
    return True

# Initialize grid
original = [
    ['f','','','c','a','','g'],
    ['e','b','c','a','d','g',''],
    ['','c','','','g','f','e'],
    ['c','a','','','','e','b'],
    ['a','d','g','','e','','c'],
    ['d','g','','e','','c',''],
    ['','','','b','','a','']
]

grid = [row[:] for row in original]  # Make a copy

# Solve
if solve(grid, original):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")