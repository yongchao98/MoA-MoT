def print_grid(grid):
    return '\n'.join(','.join(row) for row in grid)

def is_valid_complete(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check minor diagonal
    for i in range(7):
        if grid[i][6-i] != 'f':
            return False
    
    return True

def is_valid_partial(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6 and letter != 'f':
        return False
    
    return True

def solve():
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
        elif grid[i][6-i] != 'f':
            return "No solution exists"
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve_recursive():
        pos = find_empty()
        if not pos:
            return is_valid_complete(grid)
        
        row, col = pos
        for letter in 'abcdefg':
            if is_valid_partial(grid, row, col, letter):
                grid[row][col] = letter
                if solve_recursive():
                    return True
                grid[row][col] = ''
        return False
    
    if solve_recursive():
        return print_grid(grid)
    return "No solution exists"

print(solve())