def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_move(grid, row, col, letter, diagonal_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal constraint
    if row + col == 6 and letter != diagonal_letter:
        return False
            
    return True

def verify_initial_grid(grid):
    # Check for conflicts in pre-filled positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                # Temporarily remove the letter and check if it's valid
                letter = grid[i][j]
                grid[i][j] = ''
                if not is_valid_move(grid, i, j, letter, 'c'):
                    return False
                grid[i][j] = letter
    return True

def solve(grid):
    if not verify_initial_grid(grid):
        return False
    
    def solve_cell(row, col):
        if row == 7:
            return True
            
        next_row = row + 1 if col == 6 else row
        next_col = 0 if col == 6 else col + 1
        
        # If cell is already filled, verify it's valid and move on
        if grid[row][col] != '':
            if row + col == 6 and grid[row][col] != 'c':
                return False
            return solve_cell(next_row, next_col)
        
        # If on diagonal, must use 'c'
        if row + col == 6:
            if is_valid_move(grid, row, col, 'c', 'c'):
                grid[row][col] = 'c'
                if solve_cell(next_row, next_col):
                    return True
                grid[row][col] = ''
            return False
        
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter, 'c'):
                grid[row][col] = letter
                if solve_cell(next_row, next_col):
                    return True
                grid[row][col] = ''
        return False
    
    return solve_cell(0, 0)

# Initialize grid
initial_grid = [
    ['g', 'b', '', '', '', 'a', 'c'],
    ['b', 'e', 'f', 'd', 'a', 'c', 'g'],
    ['', '', '', '', 'c', 'g', ''],
    ['f', '', '', '', '', '', 'e'],
    ['', 'a', 'c', '', 'b', 'e', 'f'],
    ['a', '', 'g', 'b', 'e', 'f', ''],
    ['c', 'g', 'b', '', '', 'd', '']
]

# Create a working copy
grid = [row[:] for row in initial_grid]

# Solve and print
if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")