def solve_puzzle():
    # Initial configuration
    initial = [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]
    
    # Create working grid
    grid = [[''] * 7 for _ in range(7)]
    
    def is_valid(row, col, letter):
        # Check initial configuration
        if initial[row][col] != '' and initial[row][col] != letter:
            return False
            
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
                j = 6 - i
                if grid[i][j] != '' and grid[i][j] != letter:
                    return False
                    
        return True
    
    def find_empty():
        # First try to fill cells that are in initial configuration
        for i in range(7):
            for j in range(7):
                if initial[i][j] != '' and grid[i][j] == '':
                    return (i, j)
        # Then fill remaining empty cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve():
        pos = find_empty()
        if not pos:
            return True
            
        row, col = pos
        
        # If this is a pre-filled position, try only that letter
        if initial[row][col] != '':
            letter = initial[row][col]
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
            return False
        
        # If on minor diagonal, try to match existing diagonal letter
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    letter = grid[i][j]
                    if is_valid(row, col, letter):
                        grid[row][col] = letter
                        if solve():
                            return True
                        grid[row][col] = ''
                    return False
        
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
                
        return False
    
    # First, copy all initial values
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]
    
    if solve():
        for row in grid:
            print(','.join(row))
        return True
    else:
        print("No solution exists")
        return False

solve_puzzle()