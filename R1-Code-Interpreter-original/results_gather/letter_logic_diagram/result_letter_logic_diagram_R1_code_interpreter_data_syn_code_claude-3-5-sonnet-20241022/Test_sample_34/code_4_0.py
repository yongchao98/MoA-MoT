def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_prefilled(grid, original):
    for i in range(7):
        for j in range(7):
            if original[i][j] != '' and grid[i][j] != original[i][j]:
                return False
    return True

def get_diagonal_letter(original):
    # Get the first non-empty letter on the minor diagonal
    for i in range(7):
        if original[i][6-i] != '':
            return original[i][6-i]
    return None

def is_valid_complete(grid, diag_letter):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    for i in range(7):
        if grid[i][6-i] != diag_letter:
            return False
    
    return True

def solve_puzzle():
    original = [
        ['f','','','c','a','','g'],
        ['e','b','c','a','d','g',''],
        ['','c','','','g','f','e'],
        ['c','a','','','','e','b'],
        ['a','d','g','','e','','c'],
        ['d','g','','e','','c',''],
        ['','','','b','','a','']
    ]
    
    diag_letter = get_diagonal_letter(original)
    if not diag_letter:
        return None
    
    # First, create a grid with all minor diagonal cells filled
    grid = [[''] * 7 for _ in range(7)]
    
    # Copy pre-filled cells and fill diagonal
    for i in range(7):
        for j in range(7):
            if original[i][j] != '':
                grid[i][j] = original[i][j]
            elif i + j == 6:  # minor diagonal
                grid[i][j] = diag_letter
    
    def can_place(row, col, letter):
        # Check row
        for j in range(7):
            if grid[row][j] == letter:
                return False
        
        # Check column
        for i in range(7):
            if grid[i][col] == letter:
                return False
        
        # Check if it's a diagonal position
        if row + col == 6 and letter != diag_letter:
            return False
            
        return True
    
    def solve_recursive(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return is_valid_complete(grid, diag_letter)
        
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells and diagonal cells
        if grid[row][col] != '':
            return solve_recursive(pos + 1)
        
        for letter in 'abcdefg':
            if can_place(row, col, letter):
                grid[row][col] = letter
                if solve_recursive(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve_recursive():
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    print_grid(solution)
else:
    print("No solution found")