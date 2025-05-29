def can_use_on_diagonal(grid, letter):
    # Check if 'letter' can be used on all diagonal positions
    for i in range(7):
        j = 6-i  # diagonal position
        
        # If position is filled, it must match
        if grid[i][j] != '' and grid[i][j] != letter:
            return False
            
        # Check row conflicts
        for k in range(7):
            if k != j and grid[i][k] == letter:
                return False
                
        # Check column conflicts
        for k in range(7):
            if k != i and grid[k][j] == letter:
                return False
    
    return True

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = [diagonal_letter] if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
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

# Try each possible letter for the diagonal
for diagonal_letter in 'abcdefg':
    if can_use_on_diagonal(grid, diagonal_letter):
        # Create a fresh grid
        test_grid = [row[:] for row in grid]
        
        # Fill empty diagonal positions
        for i in range(7):
            if test_grid[i][6-i] == '':
                test_grid[i][6-i] = diagonal_letter
        
        if solve(test_grid, diagonal_letter):
            print("<<<")
            for row in test_grid:
                print(','.join(row))
            print(">>>")
            exit()

print("No solution exists")