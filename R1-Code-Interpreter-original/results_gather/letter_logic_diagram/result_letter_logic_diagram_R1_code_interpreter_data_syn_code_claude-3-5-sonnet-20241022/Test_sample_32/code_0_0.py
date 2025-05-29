def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this position is on it
    if row + col == 6:
        for i in range(7):
            if i != row and grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

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
    
    def backtrack():
        empty = find_empty(grid)
        if not empty:
            return True
        
        row, col = empty
        letters = 'abcdefg'
        
        # If position is on minor diagonal, we must use the same letter as others
        if row + col == 6:
            for i in range(7):
                if grid[i][6-i] != '':
                    letters = grid[i][6-i]
                    break
            if letters == 'abcdefg':  # No letter found on diagonal yet
                letters = 'b'  # We can deduce 'b' must be on diagonal from constraints
        
        for letter in letters:
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