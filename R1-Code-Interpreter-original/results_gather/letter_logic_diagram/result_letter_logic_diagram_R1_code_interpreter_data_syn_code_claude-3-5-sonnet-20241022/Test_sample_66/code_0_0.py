def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal consistency if this cell is on minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
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
            
            # If this is on minor diagonal, ensure all empty minor diagonal cells get same letter
            if row + col == 6:
                valid = True
                for i in range(7):
                    if grid[i][6-i] == '':
                        if is_valid(grid, i, 6-i, letter):
                            grid[i][6-i] = letter
                        else:
                            valid = False
                            break
                if not valid:
                    grid[row][col] = ''
                    continue
            
            if solve_grid(grid, row, col + 1):
                return True
            
            # Backtrack
            grid[row][col] = ''
            if row + col == 6:
                for i in range(7):
                    if grid[i][6-i] == letter:
                        grid[i][6-i] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'c', '', 'e', '', '', ''],
    ['c', '', 'e', '', '', '', 'd'],
    ['a', '', 'b', 'f', '', '', 'c'],
    ['', '', '', '', '', 'c', 'a'],
    ['b', '', 'g', 'd', '', '', ''],
    ['', 'g', '', '', '', 'e', 'b'],
    ['', '', 'c', '', 'e', '', '']
]

# Solve the puzzle
if solve_grid(initial_grid):
    # Format the output
    result = []
    for row in initial_grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")