def solve():
    # Initialize grid with givens
    grid = [
        ['', 'a', 'c', '', '', 'e', 'd'],
        ['a', '', 'f', 'g', '', 'd', 'b'],
        ['c', 'f', '', '', '', '', ''],
        ['', 'g', 'e', 'd', '', '', ''],
        ['', 'e', 'd', 'b', '', '', ''],
        ['', '', '', '', '', '', ''],
        ['d', '', 'a', '', '', '', '']
    ]
    
    # First, forcefully fill minor diagonal with 'd'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'd'
        elif grid[i][6-i] != 'd':
            return None  # Invalid initial state
    
    def is_valid(row, col, val):
        if row + col == 6:  # Minor diagonal
            return val == 'd'
        
        # Check row
        for j in range(7):
            if grid[row][j] == val:
                return False
        
        # Check column
        for i in range(7):
            if grid[i][col] == val:
                return False
        
        return True
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve_recursive():
        pos = find_empty()
        if not pos:
            return True
        
        row, col = pos
        for val in 'abcdefg':
            if is_valid(row, col, val):
                grid[row][col] = val
                if solve_recursive():
                    return True
                grid[row][col] = ''
        return False
    
    if solve_recursive():
        return grid
    return None

result = solve()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")