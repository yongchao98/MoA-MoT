def solve_grid():
    # Initialize grid
    grid = [
        ['b', '', '', '', '', 'a', 'f'],
        ['', '', '', '', '', '', 'b'],
        ['', 'e', 'd', 'a', '', '', ''],
        ['e', '', '', 'f', 'b', 'g', ''],
        ['', 'a', '', 'b', 'g', '', 'e'],
        ['', '', 'b', 'g', '', '', 'd'],
        ['', '', '', '', '', 'd', '']
    ]
    
    # First, ensure all minor diagonal positions can be 'f'
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != 'f':
            return "No solution exists"
        grid[i][j] = 'f'
    
    def is_valid(row, col, num):
        # Skip check for minor diagonal as we've already filled it
        if row + col == 6:
            return num == 'f'
        
        # Check row
        for j in range(7):
            if grid[row][j] == num:
                return False
        
        # Check column
        for i in range(7):
            if grid[i][col] == num:
                return False
        
        return True
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    def solve():
        find = find_empty()
        if not find:
            return True
        
        row, col = find
        for num in 'abcdefg':
            if is_valid(row, col, num):
                grid[row][col] = num
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    if solve():
        return '\n'.join(','.join(row) for row in grid)
    return "No solution exists"

print(solve_grid())