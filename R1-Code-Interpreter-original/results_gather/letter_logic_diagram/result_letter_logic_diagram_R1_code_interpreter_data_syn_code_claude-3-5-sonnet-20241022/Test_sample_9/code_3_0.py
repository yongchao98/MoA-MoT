def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    return True

def solve_grid():
    # Create empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Fill in the given values
    initial_values = [
        ['', 'c', '', '', 'b', '', ''],
        ['c', 'f', '', '', 'a', '', ''],
        ['', '', 'b', 'a', 'e', '', 'c'],
        ['', 'b', '', '', '', '', 'f'],
        ['', '', 'e', 'g', 'c', 'f', ''],
        ['', '', 'g', '', '', 'd', ''],
        ['e', '', '', '', '', 'b', '']
    ]
    
    for i in range(7):
        for j in range(7):
            if initial_values[i][j]:
                grid[i][j] = initial_values[i][j]
    
    # Try each possible letter for the minor diagonal
    for diag_letter in 'abcdefg':
        # Create a copy of the grid
        test_grid = [row[:] for row in grid]
        
        # Check if diagonal letter is possible
        valid = True
        for i in range(7):
            if test_grid[i][6-i] and test_grid[i][6-i] != diag_letter:
                valid = False
                break
            # Check if letter already exists in row or column
            if test_grid[i][6-i] == '':
                for j in range(7):
                    if (test_grid[i][j] == diag_letter or 
                        test_grid[j][6-i] == diag_letter):
                        valid = False
                        break
        
        if not valid:
            continue
            
        # Fill diagonal with the letter
        for i in range(7):
            test_grid[i][6-i] = diag_letter
            
        # Try to solve the rest
        if solve_remaining(test_grid, 0, 0):
            return test_grid
    
    return None

def solve_remaining(grid, row, col):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
        
    # Skip if cell is already filled or is on minor diagonal
    if grid[row][col] != '' or row + col == 6:
        return solve_remaining(grid, row, col + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

solution = solve_grid()
if solution:
    print_grid(solution)
else:
    print("No solution exists")