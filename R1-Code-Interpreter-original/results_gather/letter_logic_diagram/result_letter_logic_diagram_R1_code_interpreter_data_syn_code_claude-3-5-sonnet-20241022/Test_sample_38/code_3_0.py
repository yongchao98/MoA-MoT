def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_against_original(grid, original):
    for i in range(7):
        for j in range(7):
            if original[i][j] != '' and grid[i][j] != original[i][j]:
                return False
    return True

def is_valid_state(grid):
    # Check rows
    for row in grid:
        if '' in row:
            continue
        if len(set(row)) != 7:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if '' in col:
            continue
        if len(set(col)) != 7:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    non_empty = [x for x in diag if x != '']
    if len(non_empty) > 1 and len(set(non_empty)) > 1:
        return False
    
    return True

def solve_puzzle():
    original = [
        ['a','e','','f','g','c',''],
        ['','b','','','','d','a'],
        ['b','f','g','c','d','a','e'],
        ['f','','c','','a','','b'],
        ['g','c','d','a','e','b','f'],
        ['c','','','e','b','','g'],
        ['d','','','','f','g','c']
    ]
    
    grid = [row[:] for row in original]
    
    def find_empty():
        # First try minor diagonal
        for i in range(7):
            if i + (6-i) == 6 and grid[i][6-i] == '':
                return i, 6-i
        # Then rest of grid
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    def is_safe(row, col, letter):
        # Check original constraints
        if original[row][col] != '' and original[row][col] != letter:
            return False
        
        # Check row
        if letter in grid[row]:
            return False
            
        # Check column
        if letter in [grid[i][col] for i in range(7)]:
            return False
            
        # Check minor diagonal
        if row + col == 6:
            diag = [grid[i][6-i] for i in range(7)]
            non_empty = [x for x in diag if x != '']
            if non_empty and non_empty[0] != letter:
                return False
                
        return True
    
    def backtrack():
        pos = find_empty()
        if not pos:
            return verify_against_original(grid, original) and is_valid_state(grid)
            
        row, col = pos
        
        # If on minor diagonal, try only the existing diagonal letter if any
        if row + col == 6:
            diag = [grid[i][6-i] for i in range(7)]
            non_empty = [x for x in diag if x != '']
            if non_empty:
                if is_safe(row, col, non_empty[0]):
                    grid[row][col] = non_empty[0]
                    if backtrack():
                        return True
                    grid[row][col] = ''
                return False
        
        for letter in 'abcdefg':
            if is_safe(row, col, letter):
                grid[row][col] = letter
                if backtrack():
                    return True
                grid[row][col] = ''
        
        return False
    
    if backtrack():
        print_grid(grid)
        return True
    return False

solve_puzzle()