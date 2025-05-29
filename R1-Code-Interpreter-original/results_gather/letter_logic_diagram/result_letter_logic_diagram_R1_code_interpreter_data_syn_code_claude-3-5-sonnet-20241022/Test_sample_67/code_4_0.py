def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def get_diagonal_letter(grid):
    # Check existing diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
    return None

def can_use_diagonal_letter(grid, letter):
    # Check if letter can be used in all diagonal positions
    for i in range(7):
        j = 6 - i
        # If position is filled and not with this letter, invalid
        if grid[i][j] != '' and grid[i][j] != letter:
            return False
        # Check if letter exists in same row or column
        for k in range(7):
            if k != j and grid[i][k] == letter:
                return False
            if k != i and grid[k][j] == letter:
                return False
    return True

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal constraint
    if row + col == 6:  # If on diagonal
        diagonal_letter = get_diagonal_letter(grid)
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve_grid(grid):
    # First, handle the diagonal
    diagonal_letter = get_diagonal_letter(grid)
    if not diagonal_letter:
        # Try each letter for diagonal
        for letter in 'abcdefg':
            if can_use_diagonal_letter(grid, letter):
                # Fill diagonal with this letter
                for i in range(7):
                    j = 6 - i
                    if grid[i][j] == '':
                        grid[i][j] = letter
                diagonal_letter = letter
                break
    
    if not diagonal_letter:
        return False
    
    # Now solve the rest
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid with deep copy
grid = [
    ['', 'c', 'b', '', '', 'd', ''],
    ['c', '', 'g', 'e', 'd', 'f', 'a'],
    ['b', '', '', '', '', '', 'c'],
    ['', '', '', '', 'a', '', ''],
    ['', '', 'f', '', '', '', ''],
    ['', 'f', '', 'c', '', 'g', ''],
    ['f', '', '', 'b', 'g', '', '']
]

# Create a fresh copy of the grid
grid = [[cell for cell in row] for row in grid]

# Try each possible diagonal letter
for diagonal_letter in 'abcdefg':
    test_grid = [[cell for cell in row] for row in grid]
    if can_use_diagonal_letter(test_grid, diagonal_letter):
        # Fill diagonal
        for i in range(7):
            j = 6 - i
            if test_grid[i][j] == '':
                test_grid[i][j] = diagonal_letter
        
        if solve_grid(test_grid):
            result = []
            for row in test_grid:
                result.append(','.join(row))
            print('<<<')
            for row in result:
                print(row)
            print('>>>')
            break
else:
    print("No solution exists")