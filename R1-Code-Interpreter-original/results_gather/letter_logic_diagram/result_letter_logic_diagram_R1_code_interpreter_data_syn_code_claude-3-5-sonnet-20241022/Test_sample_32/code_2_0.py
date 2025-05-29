def solve_grid():
    # Initialize grid with given values
    grid = [
        ['', 'b', '', 'c', 'd', 'f', 'a'],
        ['', '', 'c', '', '', '', 'g'],
        ['', 'c', '', '', '', 'g', ''],
        ['c', 'd', '', 'a', '', '', 'e'],
        ['d', '', 'a', '', 'b', 'e', ''],
        ['f', 'a', '', '', '', '', ''],
        ['a', '', 'b', 'e', 'c', '', 'f']
    ]
    
    # Fill all diagonal positions with 'a'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'a'
    
    def is_valid(grid, row, col, letter):
        # Skip check if it's a pre-filled position
        if grid[row][col] != '':
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

    def find_empty(grid):
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None

    def backtrack():
        empty = find_empty(grid)
        if not empty:
            return True
        
        row, col = empty
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack():
                    return True
                grid[row][col] = ''
        return False
    
    if backtrack():
        result = []
        for row in grid:
            result.append(','.join(row))
        print('\n'.join(result))
    else:
        print("No solution exists")

solve_grid()