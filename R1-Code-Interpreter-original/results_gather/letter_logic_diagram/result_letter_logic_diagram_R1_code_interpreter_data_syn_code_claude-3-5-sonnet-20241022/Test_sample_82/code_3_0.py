def is_valid(grid, row, col, letter, diagonal_letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check diagonal constraint
    if row + col == 6 and letter != diagonal_letter:
        return False
        
    return True

def verify_solution(grid, diagonal_letter):
    # Verify rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify diagonal
    for i in range(7):
        if grid[i][6-i] != diagonal_letter:
            return False
    
    return True

def find_empty(grid):
    # First check diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then check other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return verify_solution(grid, diagonal_letter)
    
    row, col = empty
    letters = [diagonal_letter] if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

def can_use_diagonal(initial_grid, diagonal_letter):
    # Check if diagonal letter conflicts with pre-filled cells
    for i in range(7):
        if initial_grid[i][6-i] != '' and initial_grid[i][6-i] != diagonal_letter:
            return False
    return True

# Initial grid
initial_grid = [
    ['f','','e','','','',''],
    ['','e','','d','','','f'],
    ['','','','','g','f','a'],
    ['','d','b','','','a',''],
    ['d','b','','','a','e',''],
    ['b','g','','a','e','c','d'],
    ['','f','a','','','d','']
]

# Try each possible diagonal letter
found_solution = False
for diagonal_letter in 'abcdefg':
    if can_use_diagonal(initial_grid, diagonal_letter):
        grid = [row[:] for row in initial_grid]
        # Fill all diagonal positions with the diagonal letter
        for i in range(7):
            if grid[i][6-i] == '':
                grid[i][6-i] = diagonal_letter
        if solve_grid(grid, diagonal_letter):
            found_solution = True
            for row in grid:
                print(','.join(row))
            break

if not found_solution:
    print("No solution exists")