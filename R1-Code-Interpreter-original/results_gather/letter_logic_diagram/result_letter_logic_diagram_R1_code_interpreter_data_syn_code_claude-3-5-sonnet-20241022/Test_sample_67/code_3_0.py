def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def get_diagonal_positions():
    return [(i, 6-i) for i in range(7)]

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Special handling for diagonal positions
    diagonal_positions = get_diagonal_positions()
    if (row, col) in diagonal_positions:
        for pos in diagonal_positions:
            if grid[pos[0]][pos[1]] != '' and grid[pos[0]][pos[1]] != letter:
                return False
    
    return True

def solve_grid(grid):
    # First, determine what letter should be on the diagonal
    diagonal_positions = get_diagonal_positions()
    diagonal_letter = None
    
    # Check if any diagonal position is already filled
    for pos in diagonal_positions:
        if grid[pos[0]][pos[1]] != '':
            diagonal_letter = grid[pos[0]][pos[1]]
            break
    
    # If no diagonal letter is set, try each possibility
    if not diagonal_letter:
        for letter in 'abcdefg':
            valid = True
            # Check if this letter can be used for all diagonal positions
            for pos in diagonal_positions:
                row, col = pos
                temp_grid = [row[:] for row in grid]
                temp_grid[row][col] = letter
                if not is_valid(temp_grid, row, col, letter):
                    valid = False
                    break
            if valid:
                diagonal_letter = letter
                break
    
    if diagonal_letter is None:
        return False
    
    # Fill all diagonal positions first
    for pos in diagonal_positions:
        if grid[pos[0]][pos[1]] == '':
            grid[pos[0]][pos[1]] = diagonal_letter
    
    # Then solve the rest
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', 'c', 'b', '', '', 'd', ''],
    ['c', '', 'g', 'e', 'd', 'f', 'a'],
    ['b', '', '', '', '', '', 'c'],
    ['', '', '', '', 'a', '', ''],
    ['', '', 'f', '', '', '', ''],
    ['', 'f', '', 'c', '', 'g', ''],
    ['f', '', '', 'b', 'g', '', '']
]

# Create a deep copy of the grid for each row
grid = [[cell for cell in row] for row in grid]

if solve_grid(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")