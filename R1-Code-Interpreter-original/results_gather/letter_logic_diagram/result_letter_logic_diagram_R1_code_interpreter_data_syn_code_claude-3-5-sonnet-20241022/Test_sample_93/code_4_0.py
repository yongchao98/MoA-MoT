def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    return True

def solve(grid):
    # First, identify the letter that must be on the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letter = grid[i][6-i]
            break
    
    if not diagonal_letter:
        return False
        
    # Force all diagonal positions to have this letter
    for i in range(7):
        if grid[i][6-i] == '':
            if not is_valid_move(grid, i, 6-i, diagonal_letter):
                return False
            grid[i][6-i] = diagonal_letter
        elif grid[i][6-i] != diagonal_letter:
            return False
    
    def solve_remaining(row, col):
        if row == 7:
            return True
            
        next_row = row + 1 if col == 6 else row
        next_col = 0 if col == 6 else col + 1
        
        # Skip if cell is already filled
        if grid[row][col] != '':
            return solve_remaining(next_row, next_col)
            
        # Skip if this is a diagonal position (already filled)
        if row + col == 6:
            return solve_remaining(next_row, next_col)
        
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter):
                grid[row][col] = letter
                if solve_remaining(next_row, next_col):
                    return True
                grid[row][col] = ''
        return False
    
    return solve_remaining(0, 0)

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