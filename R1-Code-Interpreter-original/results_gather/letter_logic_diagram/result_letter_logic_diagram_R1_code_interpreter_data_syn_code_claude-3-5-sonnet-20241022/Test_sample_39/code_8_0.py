def print_grid(grid):
    for row in grid:
        print(','.join(str(x) for x in row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve_with_fixed_cells():
    # Initial configuration with fixed cells
    fixed = [
        ['', 'b', '', '', '', '', ''],
        ['b', 'a', 'g', 'd', '', '', 'f'],
        ['', '', 'd', '', '', '', 'b'],
        ['g', '', '', 'e', '', '', ''],
        ['', '', '', '', 'b', '', ''],
        ['c', '', 'f', '', 'a', '', 'd'],
        ['e', 'f', '', 'a', 'g', '', '']
    ]
    
    # Create working grid with fixed values
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # First, verify and place all fixed values
    for i in range(7):
        for j in range(7):
            if fixed[i][j] != '':
                if not is_valid(grid, i, j, fixed[i][j]):
                    return None
                grid[i][j] = fixed[i][j]
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and fixed[i][j] == '':
                    return (i, j)
        return None
    
    def solve():
        pos = find_empty()
        if not pos:
            return True
        
        row, col = pos
        
        # If on minor diagonal, check for existing diagonal value
        diagonal_value = None
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    diagonal_value = grid[i][j]
                    break
        
        # Try each possible letter
        possible_letters = [diagonal_value] if diagonal_value else 'abcdefg'
        for letter in possible_letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve():
        return grid
    return None

# Solve and print result
print('<<<')
result = solve_with_fixed_cells()
if result:
    print_grid(result)
else:
    print("No solution exists")
print('>>>')