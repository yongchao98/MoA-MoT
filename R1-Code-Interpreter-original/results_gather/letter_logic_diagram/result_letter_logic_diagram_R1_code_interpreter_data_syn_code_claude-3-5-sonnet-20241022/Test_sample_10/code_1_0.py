def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_minor_diagonal_letter(grid):
    # Check pre-filled positions on minor diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    
    # If we have exactly one letter, that must be our diagonal letter
    if len(diagonal_letters) == 1:
        return list(diagonal_letters)[0]
    
    # If we have 'f' in a critical position, use it
    if 'f' in diagonal_letters:
        return 'f'
    
    return None

def is_valid_placement(grid, row, col, letter, diagonal_letter):
    # If this is a minor diagonal position, must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
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

def find_empty_cell(grid):
    # First try to fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, diagonal_letter):
    # Find empty cell
    empty = find_empty_cell(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If this is a minor diagonal position, only try diagonal_letter
    if row + col == 6:
        if is_valid_placement(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    else:
        # Try each letter except diagonal_letter (unless necessary)
        for letter in 'abcdefg':
            if is_valid_placement(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if solve(grid, diagonal_letter):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['d','','','e','b','','g'],
    ['','','e','b','','','d'],
    ['','','b','','','','c'],
    ['','','a','','d','c','f'],
    ['','a','g','d','','','e'],
    ['','','','','','',''],
    ['','','','f','','','a']
]

# Try each possible diagonal letter
for diagonal_letter in 'abcdefg':
    # Create a fresh copy of the grid
    grid = [row[:] for row in initial_grid]
    
    # Check if this diagonal letter conflicts with existing positions
    valid = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diagonal_letter:
            valid = False
            break
    
    if not valid:
        continue
    
    # Fill all empty diagonal positions with the chosen letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    # Try to solve with this diagonal letter
    if solve(grid, diagonal_letter):
        # Verify solution
        valid = True
        # Check rows and columns
        for i in range(7):
            row_letters = set(grid[i])
            col_letters = set(grid[j][i] for j in range(7))
            if len(row_letters) != 7 or len(col_letters) != 7:
                valid = False
                break
        
        if valid:
            result = []
            for row in grid:
                result.append(','.join(row))
            print('\n'.join(result))
            break