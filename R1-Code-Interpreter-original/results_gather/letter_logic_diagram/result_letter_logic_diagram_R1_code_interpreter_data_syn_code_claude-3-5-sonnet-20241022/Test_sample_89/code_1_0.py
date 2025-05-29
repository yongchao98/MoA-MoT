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

def check_initial_match(grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def is_valid_placement(grid, row, col, letter):
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

# Start with copying all pre-filled values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            if not is_valid_placement(grid, i, j, initial_grid[i][j]):
                print("Initial configuration is invalid!")
                exit()
            grid[i][j] = initial_grid[i][j]

# Find the letter that should be on the minor diagonal
diagonal_letter = None
for i in range(7):
    j = 6 - i
    if initial_grid[i][j] != '':
        diagonal_letter = initial_grid[i][j]
        break

if diagonal_letter is None:
    print("No diagonal letter found in initial configuration!")
    exit()

# Fill the minor diagonal
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = diagonal_letter

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return check_initial_match(grid)
        
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
        
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
            
    return False

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")