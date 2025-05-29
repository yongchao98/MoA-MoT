def verify_initial_state(grid):
    # First verify all pre-filled positions
    initial = [
        ['', '', '', 'e', 'd', 'c', ''],
        ['f', '', 'e', '', '', '', ''],
        ['', '', 'd', '', 'g', '', ''],
        ['e', '', 'c', '', 'b', 'f', ''],
        ['d', '', '', 'b', '', 'a', 'e'],
        ['c', 'g', 'b', '', 'a', 'e', 'd'],
        ['', 'b', 'f', '', '', 'd', '']
    ]
    
    # Create new grid with g on diagonal
    new_grid = [[''] * 7 for _ in range(7)]
    
    # Fill diagonal with g
    for i in range(7):
        new_grid[i][6-i] = 'g'
    
    # Copy pre-filled values that aren't on diagonal
    for i in range(7):
        for j in range(7):
            if i + j != 6 and initial[i][j] != '':
                if new_grid[i][j] != '' and new_grid[i][j] != initial[i][j]:
                    return None  # Conflict found
                new_grid[i][j] = initial[i][j]
    
    return new_grid

def is_valid(grid, row, col, letter):
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
            if i + j != 6 and grid[i][j] == '':  # Skip diagonal positions
                return (i, j)
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for letter in 'abcdefg':
        if letter != 'g' and is_valid(grid, row, col, letter):  # Don't use g except on diagonal
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Create initial grid with g on diagonal
grid = verify_initial_state([[''] * 7 for _ in range(7)])

if grid and solve(grid):
    print("<<<")
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print(result.strip())
    print(">>>")
else:
    print("No valid solution exists")