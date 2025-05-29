def print_solution(grid):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")

initial_grid = [
    ['', '', 'a', 'g', '', 'b', 'c'],
    ['', 'a', '', '', '', 'c', 'd'],
    ['a', 'g', 'e', 'b', '', '', ''],
    ['g', '', 'b', '', 'd', 'f', 'a'],
    ['e', 'b', '', '', 'f', '', 'g'],
    ['', '', 'd', 'f', 'a', '', ''],
    ['', 'd', 'f', 'a', 'g', 'e', 'b']
]

def is_valid(grid, row, col, letter):
    # Check if this would violate pre-filled cells
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
        
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
            
    return True

def find_empty(grid):
    # First fill pre-filled cells
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] == '':
                return (i, j, True)  # True indicates this is a pre-filled cell
                
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
        
    row, col, is_prefilled = empty
    
    if is_prefilled:
        # For pre-filled cells, only try the specified letter
        letter = initial_grid[row][col]
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    else:
        # For other cells, try all possible letters
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
                
    return False

# Initialize empty grid
grid = [[''] * 7 for _ in range(7)]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")