def solve_grid():
    # Initialize grid
    grid = [
        ['b', '', '', '', '', 'a', 'f'],
        ['', '', '', '', '', 'f', 'b'],  # Force f on diagonal
        ['', 'e', 'd', 'a', 'f', '', ''],  # Force f on diagonal
        ['e', '', '', 'f', 'b', 'g', ''],
        ['', 'a', '', 'b', 'g', 'f', 'e'],  # Force f on diagonal
        ['', '', 'b', 'g', '', 'f', 'd'],  # Force f on diagonal
        ['f', '', '', '', '', 'd', '']  # Force f on diagonal
    ]
    
    def is_valid(row, col, letter):
        # Minor diagonal constraint
        if row + col == 6:
            return letter == 'f'
            
        # Row constraint
        for j in range(7):
            if grid[row][j] == letter:
                return False
                
        # Column constraint
        for i in range(7):
            if grid[i][col] == letter:
                return False
                
        return True
    
    def find_empty():
        # First try to fill non-diagonal cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and (i + j != 6):
                    return (i, j)
        return None
    
    def solve():
        pos = find_empty()
        if not pos:
            return True
            
        row, col = pos
        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    # First verify/fill minor diagonal with f
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != 'f':
            return "No solution exists"
        grid[i][j] = 'f'
    
    if solve():
        return '\n'.join(','.join(row) for row in grid)
    return "No solution exists"

print(solve_grid())