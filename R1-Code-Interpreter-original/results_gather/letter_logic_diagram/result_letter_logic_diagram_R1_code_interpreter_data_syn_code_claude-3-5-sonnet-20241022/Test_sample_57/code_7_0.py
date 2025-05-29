def get_possible_diagonal_letters(initial_grid):
    # Get letters that don't conflict with initial values on the minor diagonal
    used = set()
    possible = set('abcdefg')
    
    # Check initial values on minor diagonal
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '':
            return {initial_grid[i][j]}
    
    # Check what letters are already used in rows/columns of empty diagonal positions
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] == '':
            # Check row
            for k in range(7):
                if initial_grid[i][k] != '':
                    used.add(initial_grid[i][k])
            # Check column
            for k in range(7):
                if initial_grid[k][j] != '':
                    used.add(initial_grid[k][j])
    
    return possible - used

def is_valid(grid, row, col, letter, diagonal_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal
    if row + col == 6 and letter != diagonal_letter:
        return False
    
    return True

def solve_with_diagonal(grid, initial_grid, diagonal_letter, pos=0):
    if pos == 49:
        return True
    
    row = pos // 7
    col = pos % 7
    
    # If cell is pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve_with_diagonal(grid, initial_grid, diagonal_letter, pos + 1)
    
    # If on minor diagonal, must use diagonal_letter
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve_with_diagonal(grid, initial_grid, diagonal_letter, pos + 1):
                return True
        grid[row][col] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve_with_diagonal(grid, initial_grid, diagonal_letter, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Try each possible diagonal letter
possible_diagonal_letters = get_possible_diagonal_letters(initial_grid)
solution_found = False

for diagonal_letter in possible_diagonal_letters:
    grid = [[''] * 7 for _ in range(7)]
    
    # First fill the minor diagonal
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] == '':
            grid[i][j] = diagonal_letter
        else:
            grid[i][j] = initial_grid[i][j]
    
    # Try to solve with this diagonal letter
    if solve_with_diagonal(grid, initial_grid, diagonal_letter):
        solution_found = True
        for row in grid:
            print(','.join(row))
        break

if not solution_found:
    print("No solution exists")