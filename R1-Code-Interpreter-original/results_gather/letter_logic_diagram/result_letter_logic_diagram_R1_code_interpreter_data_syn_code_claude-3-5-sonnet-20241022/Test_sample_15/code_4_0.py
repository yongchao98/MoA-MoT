def get_diagonal_constraints(grid):
    # Return what letters are already used in rows and columns of diagonal positions
    used_in_rows = [set() for _ in range(7)]
    used_in_cols = [set() for _ in range(7)]
    diagonal_cells = set()
    
    # Fill in what we know
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                used_in_rows[i].add(grid[i][j])
                used_in_cols[j].add(grid[i][j])
                if i + j == 6:  # on diagonal
                    diagonal_cells.add(grid[i][j])
    
    # For each diagonal position, find possible letters
    possible_letters = set('abcdefg')
    for i in range(7):
        j = 6-i  # diagonal position
        # Letter can't be used in this row or column
        possible_letters -= used_in_rows[i]
        possible_letters -= used_in_cols[j]
    
    # Must also include any pre-filled diagonal cells
    if diagonal_cells:
        possible_letters &= diagonal_cells
    
    return possible_letters

def is_valid_row(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
            
    return True

def solve_row(grid, row, used_in_row):
    if row >= 7:
        return True
        
    # If row is complete, move to next row
    if '' not in grid[row]:
        return solve_row(grid, row + 1, set())
        
    # Find first empty cell in row
    for col in range(7):
        if grid[row][col] == '':
            break
    else:
        return solve_row(grid, row + 1, set())
    
    # Try each possible letter
    for letter in 'abcdefg':
        if letter not in used_in_row and is_valid_row(grid, row, col, letter):
            grid[row][col] = letter
            used_in_row.add(letter)
            if solve_row(grid, row, used_in_row):
                return True
            grid[row][col] = ''
            used_in_row.remove(letter)
            
    return False

# Initial grid
grid = [
    ['c','','','','a','','b'],
    ['g','e','f','a','','',''],
    ['','','','d','','','g'],
    ['f','a','','','c','g',''],
    ['','d','','c','','e','f'],
    ['d','','','','','f',''],
    ['','','','','','','']
]

# First, find what letter must be on the diagonal
possible_diagonal = get_diagonal_constraints(grid)
if len(possible_diagonal) == 0:
    print("No solution possible!")
    exit()

# Try each possible diagonal letter
for diag_letter in possible_diagonal:
    # Create a fresh grid
    test_grid = [row[:] for row in grid]
    
    # Fill diagonal with this letter
    for i in range(7):
        if test_grid[i][6-i] == '':
            test_grid[i][6-i] = diag_letter
    
    # Now solve row by row
    if solve_row(test_grid, 0, set(test_grid[0])):
        print("<<<")
        for row in test_grid:
            print(','.join(row))
        print(">>>")
        exit()

print("No solution found")