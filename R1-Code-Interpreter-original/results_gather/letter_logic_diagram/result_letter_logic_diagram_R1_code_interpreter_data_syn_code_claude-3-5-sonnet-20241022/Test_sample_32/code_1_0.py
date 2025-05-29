def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must be the same letter
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

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
    
    def find_empty(grid):
        # First fill diagonal positions
        for i in range(7):
            if grid[i][6-i] == '':
                return i, 6-i
        # Then fill other positions
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
        letters = 'abcdefg'
        
        # If this is a diagonal position, use the existing diagonal letter if any
        if row + col == 6:
            for i in range(7):
                if grid[i][6-i] != '':
                    letters = grid[i][6-i]
                    break
        
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