def solve_grid():
    # Initialize grid with given values
    grid = [
        ['b', '', '', '', '', 'a', 'f'],
        ['', '', '', '', '', '', 'b'],
        ['', 'e', 'd', 'a', '', '', ''],
        ['e', '', '', 'f', 'b', 'g', ''],
        ['', 'a', '', 'b', 'g', '', 'e'],
        ['', '', 'b', 'g', '', '', 'd'],
        ['', '', '', '', '', 'd', '']
    ]
    
    # First, fill all minor diagonal positions with 'f'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'f'
    
    def is_valid(row, col, num):
        # Check row
        if num in grid[row]:
            return False
            
        # Check column
        if num in [grid[i][col] for i in range(7)]:
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