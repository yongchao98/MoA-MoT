def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_original(grid, original):
    for i in range(7):
        for j in range(7):
            if original[i][j] != '' and grid[i][j] != original[i][j]:
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
            
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

# Original puzzle
original = [
    ['a','e','','f','g','c',''],
    ['','b','','','','d','a'],
    ['b','f','g','c','d','a','e'],
    ['f','','c','','a','','b'],
    ['g','c','d','a','e','b','f'],
    ['c','','','e','b','','g'],
    ['d','','','','f','g','c']
]

def solve():
    # Initialize empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # First, copy all pre-filled values
    for i in range(7):
        for j in range(7):
            if original[i][j] != '':
                grid[i][j] = original[i][j]
    
    def find_minor_diagonal_letter():
        for i in range(7):
            if grid[i][6-i] != '':
                return grid[i][6-i]
        return None
    
    def find_empty():
        # First fill minor diagonal
        for i in range(7):
            if i + (6-i) == 6 and grid[i][6-i] == '':
                return (i, 6-i)
        # Then other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve_recursive():
        pos = find_empty()
        if not pos:
            return verify_original(grid, original)
        
        row, col = pos
        
        # If this is a pre-filled position in original, skip it
        if original[row][col] != '':
            grid[row][col] = original[row][col]
            return solve_recursive()
        
        # If on minor diagonal, use the existing diagonal letter
        if row + col == 6:
            diag_letter = find_minor_diagonal_letter()
            if diag_letter:
                if is_valid(grid, row, col, diag_letter):
                    grid[row][col] = diag_letter
                    if solve_recursive():
                        return True
                    grid[row][col] = ''
                return False
        
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_recursive():
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve_recursive():
        print_grid(grid)
        return True
    else:
        print("No solution exists")
        return False

solve()