def print_solution(grid):
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print("<<<")
    print(result.strip())
    print(">>>")

# Initial grid with 'g' forced on the minor diagonal
grid = [
    ['', '', '', 'e', 'd', 'c', 'g'],  # Force g
    ['f', '', 'e', '', '', 'g', ''],   # Force g
    ['', '', 'd', '', 'g', '', ''],    # Already g
    ['e', '', 'c', 'g', 'b', 'f', ''], # Force g
    ['d', '', 'g', 'b', '', 'a', 'e'], # Force g
    ['c', 'g', 'b', '', 'a', 'e', 'd'], # Already g
    ['g', 'b', 'f', '', '', 'd', '']    # Force g
]

def is_valid(grid, row, col, letter):
    # Skip if trying to change a pre-filled cell
    if grid[row][col] != '':
        return False
    
    # Skip if trying to change diagonal 'g'
    if row + col == 6 and letter != 'g':
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
                return (i, j)
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
        
    row, col = pos
    letters = 'abcdefg'
    if row + col == 6:  # If on diagonal, only try 'g'
        letters = 'g'
        
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
            
    return False

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")