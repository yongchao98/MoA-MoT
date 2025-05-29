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
    
    def is_valid(row, col, letter):
        # Check row uniqueness
        for j in range(7):
            if grid[row][j] == letter:
                return False
        
        # Check column uniqueness
        for i in range(7):
            if grid[i][col] == letter:
                return False
        
        # Check minor diagonal
        if row + col == 6:
            return letter == 'f'
        
        return True
    
    def find_empty():
        # First fill minor diagonal
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i)
        # Then fill other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve():
        pos = find_empty()
        if not pos:
            # Verify solution
            for i in range(7):
                # Check row uniqueness
                if len(set(grid[i])) != 7:
                    return False
                # Check column uniqueness
                col = [grid[j][i] for j in range(7)]
                if len(set(col)) != 7:
                    return False
                # Check minor diagonal
                if grid[i][6-i] != 'f':
                    return False
            return True
        
        row, col = pos
        letters = ['f'] if row + col == 6 else 'abcdefg'
        
        for letter in letters:
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    if solve():
        return '\n'.join(','.join(row) for row in grid)
    return "No solution exists"

print(solve_grid())