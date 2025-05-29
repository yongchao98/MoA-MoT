def validate_and_solve():
    # Initial grid
    grid = [
        ['', 'a', 'c', '', '', 'e', 'd'],
        ['a', '', 'f', 'g', '', 'd', 'b'],
        ['c', 'f', '', '', '', '', ''],
        ['', 'g', 'e', 'd', '', '', ''],
        ['', 'e', 'd', 'b', '', '', ''],
        ['', '', '', '', '', '', ''],
        ['d', '', 'a', '', '', '', '']
    ]
    
    # First, verify that all given positions on minor diagonal are 'd'
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != 'd':
            return None  # Invalid initial grid
    
    def is_valid(row, col, letter):
        # Minor diagonal must be 'd'
        if row + col == 6 and letter != 'd':
            return False
            
        # Check row
        for j in range(7):
            if grid[row][j] == letter:
                return False
                
        # Check column
        for i in range(7):
            if grid[i][col] == letter:
                return False
                
        return True
    
    def find_next_empty():
        # First fill minor diagonal
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i)
        # Then other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve():
        pos = find_next_empty()
        if not pos:
            return True
            
        row, col = pos
        
        # If on minor diagonal, only try 'd'
        if row + col == 6:
            if is_valid(row, col, 'd'):
                grid[row][col] = 'd'
                if solve():
                    return True
                grid[row][col] = ''
        else:
            # Try each possible letter
            for letter in 'abcdefg':
                if is_valid(row, col, letter):
                    grid[row][col] = letter
                    if solve():
                        return True
                    grid[row][col] = ''
        return False
    
    # First fill all minor diagonal with 'd'
    for i in range(7):
        if grid[i][6-i] == '':
            if not is_valid(i, 6-i, 'd'):
                return None
            grid[i][6-i] = 'd'
    
    if solve():
        return grid
    return None

result = validate_and_solve()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")